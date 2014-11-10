#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "fdtdTM_upml.h"
#include "field.h"
#include "models.h"
#include "ntffTM.h"
#include "function.h"

//Ez(i    , j    ) -> Ez[i,j];
//Hx(i    , j+0.5) -> Hx[i,j];
//Hy(i+0.5, j    ) -> Hy[i,j];

static dcomplex *Ez = NULL;
static dcomplex *Jz = NULL;
static dcomplex *Dz = NULL;

static dcomplex *Hx = NULL;
static dcomplex *Mx = NULL;
static dcomplex *Bx = NULL;

static dcomplex *Hy = NULL;
static dcomplex *My = NULL;
static dcomplex *By = NULL;

//--------------------for NTFF --------------------//
static dcomplex *Ux,*Uy,*Wz;

static double *C_JZ = NULL, *C_MX = NULL, *C_MY = NULL;
static double *C_JZHXHY = NULL, *C_MXEZ = NULL, *C_MYEZ = NULL;
static double *C_DZ = NULL, *C_BX = NULL, *C_BY = NULL;
static double *C_DZJZ0 = NULL, *C_DZJZ1 = NULL;
static double *C_BXMX0 = NULL, *C_BXMX1 = NULL;
static double *C_BYMY0 = NULL, *C_BYMY1 = NULL;

static double *EPS_EZ = NULL, *EPS_HX = NULL, *EPS_HY = NULL;

static void freeMemories(void);
static void allocateMemories(void);
static void setCoefficient(void);

static void update(void);
static void finish(void);
static void reset(void);
static void init(void);

static void calcJD(void);
static void calcE(void);
static void calcMB(void);
static void calcH(void);

//Update
static void update(void)
{
  calcMB();  
  calcH();
  
  calcJD();
  calcE();
  
//  field_scatteredWave(Ez, EPS_EZ, 0, 0); //Ezは格子点上に配置されているので,ずれは(0,0)
  field_scatteredPulse(Ez, EPS_EZ, 0, 0, 1.0); //Ezは格子点上に配置されているので,ずれは(0,0)  

  ntffTM_TimeCalc(Hx,Hy,Ez,Ux,Uy,Wz);
}

//Initialize
static void init()
{
  allocateMemories();
  setCoefficient();

  ntffTM_init();
}

//Finish
static void finish()
{
  reset();
  freeMemories();
}

//Reset -> メモリとかは解放せず, もう一度シミュレーションを行う
static void reset()
{
  ntffTM_TimeOutput(Ux,Uy,Wz);

  /*  
      dcomplex res[360];
      ntffTM_Frequency(Hx, Hy, Ez, res);
      FILE *fp = fopen("res_5.txt", "w");
      for(int i=0; i<360; i++)
      fprintf(fp,"%.18lf\n", cnorm(res[i]));
  */

//計算領域の初期化
  memset(Hx, 0, sizeof(double complex)*N_CELL);
  memset(Hy, 0, sizeof(double complex)*N_CELL);
  memset(Ez, 0, sizeof(double complex)*N_CELL);

  memset(Mx, 0, sizeof(double complex)*N_CELL);
  memset(My, 0, sizeof(double complex)*N_CELL);
  memset(Jz, 0, sizeof(double complex)*N_CELL);

  memset(Bx, 0, sizeof(double complex)*N_CELL);
  memset(By, 0, sizeof(double complex)*N_CELL);
  memset(Dz, 0, sizeof(double complex)*N_CELL);

  int size = sizeof(dcomplex)*field_getNTFFInfo().arraySize * 360;
  memset(Ux, 0, size);
  memset(Uy, 0, size);
  memset(Wz, 0, size);
}

//:public
void (* fdtdTM_upml_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTM_upml_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTM_upml_getReset(void))(void)
{
  return reset;
}

void (* fdtdTM_upml_getInit(void))(void)
{
  return init;
}

double complex* fdtdTM_upml_getHx(void){
  return Hx;
}

double complex* fdtdTM_upml_getHy(void){
  return Hy;
}

double complex* fdtdTM_upml_getEz(void){
  return Ez;
}

double* fdtdTM_upml_getEps()
{
  return EPS_EZ;
}

//calculate J and D
static void calcJD()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      double complex nowJz = Jz[k];
      Jz[k] = C_JZ[k]*Jz[k] + C_JZHXHY[k]*(+Hy[k] - Hy[k-fInfo_s.N_PY] - Hx[k] + Hx[k-1] );
      Dz[k] = C_DZ[k]*Dz[k] + C_DZJZ1[k]*Jz[k] - C_DZJZ0[k]*nowJz;
    }
  }
}

//calculate E 
static void calcE()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      Ez[k] = Dz[k]/EPS_EZ[k];
    }
  }
}

//calculate M and B
static void calcMB()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      double complex nowMx = Mx[k];
      Mx[k] = C_MX[k]*Mx[k] - C_MXEZ[k]*(Ez[k+1] - Ez[k]);
      Bx[k] = C_BX[k]*Bx[k] + C_BXMX1[k]*Mx[k] - C_BXMX0[k]*nowMx;
    }
  }

  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      double complex nowMy = My[k];
      My[k] = C_MY[k]*My[k] - C_MYEZ[k]*(-Ez[k+fInfo_s.N_PY] + Ez[k]);
      By[k] = C_BY[k]*By[k] + C_BYMY1[k]*My[k] - C_BYMY0[k]*nowMy;
    }
  }
}

static void calcH()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      Hx[k] = Bx[k]/MU_0_S;
    }
  }
  
  for(int i=1; i<fInfo_s.N_PX-1; i++){
    for(int j=1; j<fInfo_s.N_PY-1; j++){
      const int k = field_index(i,j);
      Hy[k] = By[k]/MU_0_S;
    }
  }
}

//==================================================//
//初期化と解放
//==================================================//
static void setCoefficient(void){
  //Ez,, Hx, Hyそれぞれでσx,σyが違う(場所が違うから)
  double sig_ez_x, sig_ez_y;
  double sig_hx_x, sig_hx_y;
  double sig_hy_x, sig_hy_y;
  
  const double R = 1.0e-8;
  const double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*C_0_S/2.0/N_PML/cos(M_PI/3)*log(R); //岡田先輩の博士論文より(宇野さんの本ではcosが無い)
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  for(int i=0; i<fInfo_s.N_PX; i++){
    for(int j=0; j<fInfo_s.N_PY; j++){
      int k = field_index(i,j);
      EPS_EZ[k] = models_eps(i,j, D_XY);     
      EPS_HX[k] = models_eps(i,j+0.5, D_XY);//todo D_X, D_Yにしなくていいのか?
      EPS_HY[k] = models_eps(i+0.5,j, D_XY);
      
      sig_ez_x = sig_max*field_sigmaX(i,j);
      sig_ez_y = sig_max*field_sigmaY(i,j);

      sig_hx_x = sig_max*field_sigmaX(i,j+0.5);
      sig_hx_y = sig_max*field_sigmaY(i,j+0.5);
      
      sig_hy_x = sig_max*field_sigmaX(i+0.5,j);
      sig_hy_y = sig_max*field_sigmaY(i+0.5,j);

      double sig_z = 0; // σz is zero in
      
      //Δt = 1  Κ_i = 1
      double eps = EPSILON_0_S;
      
      C_JZ[k]     = ( 2*eps - sig_ez_x) / (2*eps + sig_ez_x);
      C_JZHXHY[k] = ( 2*eps ) / (2*eps + sig_ez_x);
      C_DZ[k]     = ( 2*eps - sig_ez_y)  / (2*eps + sig_ez_y);      
      C_DZJZ1[k]  = ( 2*eps + sig_z)/(2*eps + sig_ez_y);
      C_DZJZ0[k]  = ( 2*eps - sig_z)/(2*eps + sig_ez_y);

      C_MX[k]    = (2*eps - sig_hx_y) / (2*eps + sig_hx_y);
      C_MXEZ[k]  = (2*eps) / (2*eps + sig_hx_y);
      C_BX[k]    = (2*eps - sig_z) / (2*eps + sig_z);
      C_BXMX1[k] = (2*eps + sig_hx_x) / (2*eps + sig_z);
      C_BXMX0[k] = (2*eps - sig_hx_x) / (2*eps + sig_z);
      
      C_MY[k]    = (2*eps - sig_z) / (2*eps + sig_z);
      C_MYEZ[k]  = (2*eps) / (2*eps + sig_z);      
      C_BY[k]    = (2*eps - sig_hy_x) / (2*eps + sig_hy_x);
      C_BYMY1[k] = (2*eps + sig_hy_y) / (2*eps + sig_hy_x);
      C_BYMY0[k] = (2*eps - sig_hy_y) / (2*eps + sig_hy_x);      
    }
  }
}

static void allocateMemories()
{
  Ez = newDComplex(N_CELL);
  Dz = newDComplex(N_CELL);
  Jz = newDComplex(N_CELL);
  
  Hx = newDComplex(N_CELL);
  Mx = newDComplex(N_CELL);
  Bx = newDComplex(N_CELL);
  
  Hy = newDComplex(N_CELL);
  My = newDComplex(N_CELL);
  By = newDComplex(N_CELL);
  
  C_JZ = newDouble(N_CELL);
  C_MX = newDouble(N_CELL);
  C_MY = newDouble(N_CELL);

  C_DZ = newDouble(N_CELL);
  C_BX = newDouble(N_CELL);
  C_BY = newDouble(N_CELL);

  C_JZHXHY = newDouble(N_CELL);
  C_MXEZ = newDouble(N_CELL);
  C_MYEZ = newDouble(N_CELL);

  C_DZJZ0 = newDouble(N_CELL);
  C_DZJZ1 = newDouble(N_CELL);

  C_BXMX1 = newDouble(N_CELL);
  C_BXMX0 = newDouble(N_CELL);

  C_BYMY1 = newDouble(N_CELL);
  C_BYMY0 = newDouble(N_CELL);

  EPS_HY = newDouble(N_CELL);
  EPS_HX = newDouble(N_CELL);
  EPS_EZ = newDouble(N_CELL);

  int step = field_getNTFFInfo().arraySize;
  Ux = newDComplex(360*step);
  Uy = newDComplex(360*step);
  Wz = newDComplex(360*step);
}

static void freeMemories()
{
  delete(Hx);  delete(Hy);  delete(Ez);  
  delete(Mx);  delete(My);  delete(Jz);
  delete(Bx);  delete(By);  delete(Dz);
  delete(Ux);  delete(Uy);  delete(Wz);

  delete(C_JZ); delete(C_JZHXHY); 
  delete(C_DZ); delete(C_DZJZ1); delete(C_DZJZ0);

  delete(C_MX); delete(C_MXEZ);
  delete(C_BX); delete(C_BXMX1); delete(C_BXMX0);

  delete(C_MY); delete(C_MYEZ);
  delete(C_BY); delete(C_BYMY1); delete(C_BYMY0);
}

