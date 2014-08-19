#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include "fdtdTE_upml.h"
#include "field.h"
#include "models.h"
#include "ntffTE.h"
#include "function.h"

//Ex(i+0.5,j) -> Ex[i,j]
//Ey(i,j+0.5) -> Ey[i,j]
//Hz(i+0.5,j+0.5) -> Hz[i,j]
static double complex *Ex = NULL;
static double complex *Jx = NULL;
static double complex *Dx = NULL;

static double complex *Ey = NULL;
static double complex *Jy = NULL;
static double complex *Dy = NULL;

static double complex *Hz = NULL;
static double complex *Mz = NULL;
static double complex *Bz = NULL;

static double *C_JX = NULL, *C_JY = NULL, *C_MZ= NULL;
static double *C_JXHZ = NULL, *C_JYHZ = NULL, *C_MZEXEY= NULL;
static double *C_DX=NULL, *C_DY=NULL, *C_BZ=NULL;
static double *C_DXJX0=NULL, *C_DXJX1=NULL;
static double *C_DYJY0=NULL, *C_DYJY1=NULL;
static double *C_BZMZ0=NULL, *C_BZMZ1=NULL;

static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_HZ=NULL;

static dcomplex *Uz, *Wx, *Wy;

//------prototype--------//
static void update(void);
static void finish(void);
static void reset(void);
static void init(void);

//高速化したcalculation
static void fastCalcE(void);
static void fastCalcJD(void);
static void fastCalcH(void);
static void fastCalcMB(void);

static void calcE(void);
static void calcJD(void);
static void calcH(void);
static void calcMB(void);

static void freeMemories(void);
static void allocateMemories(void);
static void setCoefficient(void);

//:public-------------------------------//
void (* fdtdTE_upml_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTE_upml_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTE_upml_getReset(void))(void)
{
  return reset;
}

void (* fdtdTE_upml_getInit(void))(void)
{
  return init;
}

double complex* fdtdTE_upml_getEx(void){
  return Ex;
}

double complex* fdtdTE_upml_getEy(void){
  return Ey;
}

double complex* fdtdTE_upml_getHz(void)
{
  return Hz;
}

double* fdtdTE_upml_getEps()
{
  return EPS_EY;
}

//---------------------------------------------------//


//-----------------memory allocate-------------//
static void init(){  
  allocateMemories();
  setCoefficient();
  ntffTE_init();
}

//---------------------メモリの解放--------------------//
static void finish(){
  char current[512];
  getcwd(current, 512); //カレントディレクトリを保存
  char re[1024], im[1024];
  sprintf(re, "%d[deg]_Eph_r.txt", (int)field_getWaveAngle());
  sprintf(im, "%d[deg]_Eph_i.txt", (int)field_getWaveAngle());
  FILE *fpR = openFile(re);
  FILE *fpI = openFile(im);
  ntffTE_TimeOutput(Wx, Wy, Uz, fpR, fpI);
  printf("saved %s %s & %s\n", current, re, im);
  fclose(fpR);
  fclose(fpI);
  freeMemories();
}

static void reset()
{
  char re[1024], im[1024];
  sprintf(re, "%d[deg]_Eph_r.txt", (int)field_getWaveAngle());
  sprintf(im, "%d[deg]_Eph_i.txt", (int)field_getWaveAngle());
  FILE *fpR = openFile(re);
  FILE *fpI = openFile(im);
  ntffTE_TimeOutput(Wx, Wy, Uz, fpR, fpI);
  fclose(fpR);
  fclose(fpI);

  memset(Ex, 0, sizeof(double complex)*N_CELL);
  memset(Ey, 0, sizeof(double complex)*N_CELL);
  memset(Hz, 0, sizeof(double complex)*N_CELL);

  memset(Jx, 0, sizeof(double complex)*N_CELL);
  memset(Jy, 0, sizeof(double complex)*N_CELL);
  memset(Mz, 0, sizeof(double complex)*N_CELL);

  memset(Dx, 0, sizeof(double complex)*N_CELL);
  memset(Dy, 0, sizeof(double complex)*N_CELL);
  memset(Bz, 0, sizeof(double complex)*N_CELL);

  int size = sizeof(dcomplex)*field_getNTFFInfo().arraySize * 360;
  memset(Wx, 0, size);
  memset(Wy, 0, size);
  memset(Uz, 0, size);
}

/*
//Standard Scattered Wave
static void scatteredWave(double complex *p, double *eps){
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();  
  double rad = 1.0*field_getWaveAngle()*M_PI/180.0;	//ラジアン変換
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく

  //ガウシアンパルス
    double _cos = cos(rad), _sin = sin(rad);
  const double beam_width = 50;
  const double t0 = 100;
  
  int i,j;
  for(i=N_PML; i<N_X+N_PML; i++){
    for(j=N_PML; j<N_Y+N_PML; j++){
      double ikx = i*ks_cos + j*ks_sin; //k_s*(i*cos + j*sin)
      
      //p[ind(i,j)] += ray_coef*(EPSILON_0_S/eps[ind(i,j)] - 1)*( cos(ikx-w_s*time) + I*sin(ikx-w_s*time) );
      int k = field_index(i,j);
      //ガウシアンパルス
      const double r = (i*_cos+j*_sin)/C_0_S-(time-t0);
      const double gaussian_coef = exp( -pow(r/beam_width, 2 ) );
      p[k] += gaussian_coef*(EPSILON_0_S/eps[k] - 1)*cexp(I*r*w_s);    
      
    }
  }
  }*/

static inline void update(void)
{
//  fastCalcJD();
//  fastCalcE();
  calcJD();
  calcE();

  //波数から90°回転した方向に足し合わせる.
  WaveInfo_S wInfo = field_getWaveInfo_S();
  double co = cos( (wInfo.Angle_deg+90) * M_PI/ 180.0);
  double si = sin( (wInfo.Angle_deg+90) * M_PI/ 180.0);

  if(co != 0.0)
    field_scatteredPulse(Ex, EPS_EX, 0.5, 0.0, co); //Exは格子点より右に0.5ずれた位置に配置
  if(si != 0.0)
    field_scatteredPulse(Ey, EPS_EY, 0.0, 0.5, si); //Eyは格子点より上に0.5ずれた位置に配置

//  if(co != 0.0)
//  field_scatteredWave(Ex, EPS_EX, 0.5, 0.0);
//  if(si != 0.0)
//  field_scatteredWave(Ey, EPS_EY, 0.0, 0.5);
  
//  fastCalcMB();
//  fastCalcH();
  
  calcMB();
  calcH();
  
  ntffTE_TimeCalc(Ex,Ey,Hz,Wx,Wy,Uz);
}


static void fastCalcJD()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  // for(i=1..N_PX-1)
  //   for(j=1..N_PY-1) と同じ
  FAST_FOR_FOR(k,fInfo_s)
  {
    const double complex nowJx = Jx[k];
    Jx[k] = C_JX[k]*Jx[k] + C_JXHZ[k]*(Hz[k] - Hz[k-1]);
    Dx[k] = C_DX[k]*Dx[k] + C_DXJX1[k]*Jx[k] - C_DXJX0[k]*nowJx;      
  }

  FAST_FOR_FOR(k,fInfo_s)
  {
    double complex nowJy = Jy[k];
    Jy[k] = C_JY[k]*Jy[k] + C_JYHZ[k]*(-Hz[k] + Hz[k-fInfo_s.N_PY]);
    Dy[k] = C_DY[k]*Dy[k] + C_DYJY1[k]*Jy[k] - C_DYJY0[k]*nowJy;
  }

}

//高速化したcalculation(あんまり速くなってない....)
static void fastCalcE()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  FAST_FOR_FOR(k,fInfo_s)
  {
    Ex[k] = Dx[k]/EPS_EX[k];
  }

  FAST_FOR_FOR(k,fInfo_s)
  {
    Ey[k] = Dy[k]/EPS_EY[k];
  }
}
static void fastCalcMB()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  FAST_FOR_FOR(k, fInfo_s)
  {
    dcomplex nowMz = Mz[k];
    Mz[k] = C_MZ[k]*Mz[k] - C_MZEXEY[k]*(Ey[k+fInfo_s.N_PY] - Ey[k] - Ex[k+1] + Ex[k]);
    Bz[k] = C_BZ[k]*Bz[k] + C_BZMZ1[k]*Mz[k] - C_BZMZ0[k]*nowMz; 
  }
}
static void fastCalcH()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  FAST_FOR_FOR(k,fInfo_s)
  {
    Hz[k] = Bz[k]/MU_0_S;   
  }
}


//calculation
static void calcJD(void)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  for(int i=1; i<N_PX-1; i++){
    for(int j=1; j<N_PY-1; j++){
      int k = field_index(i,j);
      double complex nowJx = Jx[k];
      //ちなみに Hz[k] - Hz[k]とすると,mie散乱入射角度0のときに波が真横に進む(理由わからんしただのバグだけどおもしろいからコメント残しておく)
      Jx[k] = C_JX[k]*Jx[k] + C_JXHZ[k]*(Hz[k] - Hz[k-1]);        
      Dx[k] = C_DX[k]*Dx[k] + C_DXJX1[k]*Jx[k] - C_DXJX0[k]*nowJx;
    }
  }
  
  for(int i=1; i<N_PX-1; i++){
    for(int j=1; j<N_PY-1; j++){
      int k = field_index(i,j);
      double complex nowJy = Jy[k];
      Jy[k] = C_JY[k]*Jy[k] + C_JYHZ[k]*(-Hz[k] + Hz[k-fInfo_s.N_PY]);
      Dy[k] = C_DY[k]*Dy[k] + C_DYJY1[k]*Jy[k] - C_DYJY0[k]*nowJy;
    }
  }
}

static void calcE(void)
{
  int i,j;
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
    {
      int k = field_index(i,j);
      Ex[k] = Dx[k]/EPS_EX[k];
    }
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
    {
      int k = field_index(i,j);
      Ey[k] = Dy[k]/EPS_EY[k];
    }
}

static void calcMB(void)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=1; i<N_PX-1; i++){
    for(int j=1; j<N_PY-1; j++){
      int k = field_index(i,j);
      double complex nowMz = Mz[k];
      Mz[k] = C_MZ[k]*Mz[k] - C_MZEXEY[k]*(Ey[k+fInfo_s.N_PY] - Ey[k] - Ex[k+1] + Ex[k]);
      Bz[k] = C_BZ[k]*Bz[k] + C_BZMZ1[k]*Mz[k] - C_BZMZ0[k]*nowMz; 
    }
  }
}

static void calcH(void)
{ 
   for(int i=1; i<N_PX-1; i++)     
    for(int j=1; j<N_PY-1; j++)
    {
      int k = field_index(i,j);
      Hz[k] = Bz[k]/MU_0_S;
    }
}

//==================================================//
//初期化と解放
//==================================================//
static void allocateMemories()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  Ex = newDComplex(fInfo_s.N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Ey = newDComplex(fInfo_s.N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Hz = newDComplex(fInfo_s.N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Jx = newDComplex(fInfo_s.N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Jy = newDComplex(fInfo_s.N_CELL);  //Ey(i,j+0.5) -> Ey[i,j]
  Mz = newDComplex(fInfo_s.N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  Dx = newDComplex(fInfo_s.N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Dy = newDComplex(fInfo_s.N_CELL);  //Ey(i,j+0.5) -> Ey[i,j]
  Bz = newDComplex(fInfo_s.N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]

  C_JX = (double *)malloc(sizeof(double)*N_CELL);
  C_JY = (double *)malloc(sizeof(double)*N_CELL);
  C_MZ = (double *)malloc(sizeof(double)*N_CELL);  
  C_JXHZ = (double *)malloc(sizeof(double)*N_CELL);
  C_JYHZ = (double *)malloc(sizeof(double)*N_CELL);
  C_MZEXEY = (double *)malloc(sizeof(double)*N_CELL);

  C_DX = (double *)malloc(sizeof(double)*N_CELL);
  C_DY = (double *)malloc(sizeof(double)*N_CELL);
  C_BZ = (double *)malloc(sizeof(double)*N_CELL);  
  C_DXJX0 = (double *)malloc(sizeof(double)*N_CELL);
  C_DXJX1 = (double *)malloc(sizeof(double)*N_CELL);
  C_DYJY0 = (double *)malloc(sizeof(double)*N_CELL);
  C_DYJY1 = (double *)malloc(sizeof(double)*N_CELL);
  C_BZMZ0 = (double *)malloc(sizeof(double)*N_CELL);
  C_BZMZ1 = (double *)malloc(sizeof(double)*N_CELL);
  
  EPS_EX = (double *)malloc(sizeof(double)*N_CELL);
  EPS_EY = (double *)malloc(sizeof(double)*N_CELL);
  EPS_HZ = (double *)malloc(sizeof(double)*N_CELL);

  NTFFInfo nInfo = field_getNTFFInfo();
  Wx = newDComplex(360*nInfo.arraySize);
  Wy = newDComplex(360*nInfo.arraySize);
  Uz = newDComplex(360*nInfo.arraySize); 
}

static void setCoefficient()
{
  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_y;
  double sig_ey_x, sig_ey_y;
  double sig_hz_x, sig_hz_y;
  double R = 1.0e-8;
  double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  int i,j;
  for(i=0; i<N_PX; i++){
    for(j=0; j<N_PY; j++){
      int k = ind(i,j);
      EPS_EX[k] = models_eps(i+0.5,j, D_Y);
      EPS_EY[k] = models_eps(i,j+0.5, D_X);
      EPS_HZ[k] = 0.5*(models_eps(i+0.5,j+0.5, D_X) + models_eps(i+0.5,j+0.5, D_Y));

      sig_ex_x = sig_max*field_sigmaX(i+0.5,j);
      sig_ex_y = sig_max*field_sigmaY(i+0.5,j);
      sig_ey_x = sig_max*field_sigmaX(i,j+0.5);
      sig_ey_y = sig_max*field_sigmaY(i,j+0.5);
      sig_hz_x = sig_max*field_sigmaX(i+0.5,j+0.5);
      sig_hz_y = sig_max*field_sigmaY(i+0.5,j+0.5);

      double eps = EPSILON_0_S;
      double sig_z = 0;

      C_JX[k] = (2*eps - sig_ex_y)/(2*eps + sig_ex_y);
      C_JXHZ[k] = (2*eps)/(2*eps + sig_ex_y);
      C_DX[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_DXJX1[k] = (2*eps + sig_ex_x) / (2*eps + sig_z);
      C_DXJX0[k] = (2*eps - sig_ex_x) / (2*eps + sig_z);

      C_JY[k] = (2*eps - sig_z) / (2*eps + sig_z);
      C_JYHZ[k] = (2*eps)/(2*eps + sig_z);
      C_DY[k] = (2*eps - sig_ey_x) / (2*eps + sig_ey_x);
      C_DYJY1[k] = (2*eps + sig_ey_y) / (2*eps + sig_ey_x);
      C_DYJY0[k] = (2*eps - sig_ey_y) / (2*eps + sig_ey_x);

      C_MZ[k] = (2*eps - sig_hz_x) / (2*eps + sig_hz_x );
      C_MZEXEY[k] = (2*eps) / (2*eps + sig_hz_x);
      C_BZ[k] = (2*eps - sig_hz_y) / (2*eps + sig_hz_y);
      C_BZMZ1[k] = (2*eps + sig_z) / (2*eps + sig_hz_y);
      C_BZMZ0[k] = (2*eps - sig_z) / (2*eps + sig_hz_y);
    }
  }  
}

static void freeMemories()
{
  if(Ex != NULL){    free(Ex); Ex = NULL;}  
  if(Ey != NULL){    free(Ey); Ey = NULL;}  
  if(Hz != NULL){    free(Hz); Hz = NULL;}

  delete(Dx);
  delete(Dy);
  delete(Bz);
  
  delete(Jx);
  delete(Jy);
  delete(Mz);

  delete(Wx);
  delete(Wy);
  delete(Uz);
  
  if(C_JX!= NULL){    free(C_JX);  C_JX = NULL;}
  if(C_JXHZ!= NULL){   free(C_JXHZ); C_JXHZ = NULL;}
  if(C_DX!= NULL){   free(C_DX); C_DX = NULL;}
  if(C_DXJX0 != NULL){   free(C_DXJX0); C_DXJX0 = NULL;}
  if(C_DXJX1 != NULL){   free(C_DXJX1); C_DXJX1 = NULL;}
  
  if(C_JY!= NULL){    free(C_JY);  C_JY = NULL;}
  if(C_JYHZ!= NULL){   free(C_JYHZ); C_JYHZ = NULL;}
  if(C_DY!= NULL){   free(C_DY); C_DY = NULL;}
  if(C_DYJY0 != NULL){   free(C_DYJY0); C_DYJY0 = NULL;}
  if(C_DYJY1 != NULL){   free(C_DYJY1); C_DYJY1 = NULL;}
  
  if(C_MZ!= NULL){    free(C_MZ);  C_MZ = NULL;}
  if(C_MZEXEY!= NULL){   free(C_MZEXEY); C_MZEXEY = NULL;}
  if(C_BZ != NULL){   free(C_BZ); C_BZ = NULL;}
  if(C_BZMZ0 != NULL){   free(C_BZMZ0); C_BZMZ0 = NULL;}
  if(C_BZMZ1 != NULL){   free(C_BZMZ1); C_BZMZ1 = NULL;}
  
  if(EPS_EX != NULL)   free(EPS_EX);
  if(EPS_EY != NULL)   free(EPS_EY);
  if(EPS_HZ != NULL)   free(EPS_HZ);
}
