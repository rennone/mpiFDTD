#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "fdtdTE.h"
#include "field.h"
#include "models.h"

//   メンバ変数    //
static double complex *Ex = NULL;
static double complex *Ey = NULL;
static double complex *Hz = NULL;
static double complex *Hzx = NULL;
static double complex *Hzy = NULL;
static double *C_EX = NULL, *C_EY = NULL, *C_EXLY = NULL, *C_EYLX = NULL, *C_HZLH = NULL;
static double *C_HZX= NULL, *C_HZY= NULL, *C_HZXLX= NULL, *C_HZYLY=NULL;
static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_HZ=NULL;

//------プロトタイプ宣言--------//
static void update(void);
static void finish(void);
static void reset(void);
static void init(void);
static inline void calcE(void);
static inline void calcH(void);


//--------------public method-----------------//
//-----------------getter-----------------------//

void (* fdtdTE_getUpdate(void))(void)
{
  return update;
}

void (* fdtdTE_getFinish(void))(void)
{
  return finish;
}

void (* fdtdTE_getReset(void))(void)
{
  return reset;
}

void (* fdtdTE_getInit(void))(void)
{
  return init;
}

double complex* fdtdTE_getEx(void){
  return Ex;
}

double complex* fdtdTE_getEy(void){
  return Ey;
}

double complex* fdtdTE_getHzx(void){
  return Hzx;
}

double complex* fdtdTE_getHzy(void){
  return Hzy;
}

double complex* fdtdTE_getHz(void){
  return Hz;
}

double* fdtdTE_getEps()
{
  return EPS_EY;
}

//-------------getter--------------------//
//--------------public method-----------------//


//--------------private getter ------------------//
static inline double complex EX(const int i, const int j)
{
  return Ex[ind(i,j)];
}

static inline double complex EY(const int i, const int j)
{
  return Ey[ind(i,j)];
}

static inline double complex HZX(const int i, const int j)
{
  return Hzx[ind(i,j)];
}

static inline double complex HZY(const int i, const int j)
{
  return Hzy[ind(i,j)];
}

static inline double complex HZ(const int i, const int j)
{
  return Hz[ind(i,j)];
}

static inline double CEX(const int i, const int j)
{
  return C_EX[ind(i,j)];
}

static inline double CEY(const int i, const int j)
{
  return C_EY[ind(i,j)];
}

static inline double CEXLY(const int i, const int j)
{
  return C_EXLY[ind(i,j)];
}

static inline double CEYLX(const int i, const int j)
{
  return C_EYLX[ind(i,j)];
}

static inline double CHZLH(const int i, const int j)
{
  return C_HZLH[ind(i,j)];
}

static inline double CHZX(const int i, const int j)
{
  return C_HZX[ind(i,j)];
}

static inline double CHZY(const int i, const int j)
{
  return C_HZY[ind(i,j)];
}

static inline double CHZXLX(const int i, const int j)
{
  return C_HZXLX[ind(i,j)];
}

static inline double CHZYLY(const int i, const int j)
{
  return C_HZYLY[ind(i,j)];
}

static inline double EPSHZ(const int i, const int j)
{
  return EPS_HZ[ind(i,j)];
}

static inline double EPSEX(const int i, const int j)
{
  return EPS_EX[ind(i,j)];
}

static inline double EPSEY(const int i, const int j)
{
  return EPS_EY[ind(i,j)];
}

//---------------------------------------------------//


//-----------------領域の設定とメモリ確保-------------//
static void init()
{  
  Ex = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Ey = (double complex*)malloc(sizeof(double complex)*N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Hzy = (double complex*)malloc(sizeof(double complex)*N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]
  Hzx = (double complex*)malloc(sizeof(double complex)*N_CELL);
  Hz = (double complex*)malloc(sizeof(double complex)*N_CELL);
  
  C_EX = (double *)malloc(sizeof(double)*N_CELL);
  C_EY = (double *)malloc(sizeof(double)*N_CELL);
  C_EXLY = (double *)malloc(sizeof(double)*N_CELL);
  C_EYLX = (double *)malloc(sizeof(double)*N_CELL);
  C_HZY  = (double *)malloc(sizeof(double)*N_CELL);
  C_HZX = (double *)malloc(sizeof(double)*N_CELL);
  C_HZLH = (double *)malloc(sizeof(double)*N_CELL);
  C_HZXLX = (double *)malloc(sizeof(double)*N_CELL);
  C_HZYLY = (double *)malloc(sizeof(double)*N_CELL);

  EPS_EX = (double *)malloc(sizeof(double)*N_CELL);
  EPS_EY = (double *)malloc(sizeof(double)*N_CELL);
  EPS_HZ = (double *)malloc(sizeof(double)*N_CELL);

  memset(Ex, 0, sizeof(double complex)*N_CELL);
  memset(Ey, 0, sizeof(double complex)*N_CELL);
  memset(Hzx,0, sizeof(double complex)*N_CELL);
  memset(Hzy,0, sizeof(double complex)*N_CELL);
  memset(Hz,0, sizeof(double complex)*N_CELL);

  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_xx, sig_ex_y, sig_ex_yy;
  double sig_ey_x, sig_ey_xx, sig_ey_y, sig_ey_yy;
  double sig_hz_x, sig_hz_xx, sig_hz_y, sig_hz_yy;
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
      sig_ex_xx = MU_0_S/EPSILON_0_S * sig_ex_x;      
      sig_ex_y = sig_max*field_sigmaY(i+0.5,j);
      sig_ex_yy = MU_0_S/EPSILON_0_S * sig_ex_y;

      sig_ey_x = sig_max*field_sigmaX(i,j+0.5);
      sig_ey_xx = MU_0_S/EPSILON_0_S * sig_ey_x;      
      sig_ey_y = sig_max*field_sigmaY(i,j+0.5);
      sig_ey_yy = MU_0_S/EPSILON_0_S * sig_ey_y;

      sig_hz_x = sig_max*field_sigmaX(i+0.5,j+0.5);
      sig_hz_xx = MU_0_S/EPSILON_0_S * sig_hz_x;      
      sig_hz_y = sig_max*field_sigmaY(i+0.5,j+0.5);
      sig_hz_yy = MU_0_S/EPSILON_0_S * sig_hz_y;
      
      C_EX[k] = field_pmlCoef(EPSEX(i,j), sig_ex_y);
      C_EXLY[k] = field_pmlCoef_LXY(EPSEX(i,j), sig_ex_y);

      C_EY[k] = field_pmlCoef(EPSEY(i,j), sig_ey_x);
      C_EYLX[k] = field_pmlCoef_LXY(EPSEY(i,j), sig_ey_x);

      C_HZX[k] = field_pmlCoef(MU_0_S, sig_hz_xx);
      C_HZXLX[k] = field_pmlCoef_LXY(MU_0_S, sig_hz_xx);
      
      C_HZY[k] = field_pmlCoef(MU_0_S, sig_hz_yy);
      C_HZYLY[k] = field_pmlCoef_LXY(MU_0_S, sig_hz_yy);
    }
  }
  
}

//---------------------メモリの解放--------------------//
static void finish()
{
  reset();
  if(Ex != NULL){    free(Ex); Ex = NULL;}
  if(Ey != NULL){    free(Ey); Ey = NULL;}  
  if(Hzx != NULL){   free(Hzx); Hzx = NULL;}
  if(Hzy != NULL){   free(Hzy); Hzy = NULL;}
  
  if(C_EY!= NULL){   free(C_EY); C_EY = NULL;}
  if(C_EX!= NULL){   free(C_EX); C_EX = NULL;}
  if(C_EXLY!= NULL){   free(C_EXLY); C_EXLY = NULL;}
  if(C_EYLX!= NULL){   free(C_EYLX); C_EYLX = NULL;}
  if(C_HZLH!= NULL){   free(C_HZLH); C_HZLH = NULL;}
  if(C_HZX != NULL){   free(C_HZX); C_HZX = NULL;}
  if(C_HZY != NULL){   free(C_HZY); C_HZY = NULL;}
  if(C_HZXLX != NULL){   free(C_HZXLX); C_HZXLX = NULL;}
  if(C_HZYLY != NULL){   free(C_HZYLY); C_HZYLY = NULL;}
  
  if(EPS_EX != NULL){   free(EPS_EX); EPS_EX = NULL;}
  if(EPS_EY != NULL){   free(EPS_EY); EPS_EY = NULL;}
  if(EPS_HZ != NULL){   free(EPS_HZ); EPS_HZ = NULL;}
}

static void reset()
{
  FieldInfo fInfo = field_getFieldInfo();
  char buf[128];
  double radius_s = field_toCellUnit(1.2*fInfo.lambda_nm);
  sprintf(buf, "te_%dnm.txt",fInfo.h_u_nm);
  field_outputElliptic(buf, Ey, radius_s);

  memset(Ex, 0, sizeof(double complex)*N_CELL);
  memset(Ey, 0, sizeof(double complex)*N_CELL);
  memset(Hzx,0, sizeof(double complex)*N_CELL);
  memset(Hzy,0, sizeof(double complex)*N_CELL);
  memset(Hz,0, sizeof(double complex)*N_CELL);
}

static void update(){
  calcE();
  field_scatteredWaveNotUPML(Ey, EPS_EY, 0.0, 0.5);
  calcH();
}

//電界の計算 
static inline void calcE(void)
{
  int i,j;
  //Ex
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ex[ind(i,j)] = CEX(i,j)*EX(i,j) + CEXLY(i,j)*( HZX(i,j) - HZX(i,j-1) + HZY(i,j) - HZY(i,j-1) );
  
  //Ey
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Ey[ind(i,j)] = CEY(i,j)*EY(i,j) - CEYLX(i,j)*( HZX(i,j) - HZX(i-1,j) + HZY(i,j) - HZY(i-1,j) );
}

//磁界の計算 
static inline void calcH()
{
  int i,j;
  //Hzx
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Hzx[ind(i,j)] = CHZX(i,j)*HZX(i,j) - CHZXLX(i,j)*(EY(i+1,j)-EY(i,j) );

  //Hzy
  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Hzy[ind(i,j)] = CHZY(i,j)*HZY(i,j) + CHZYLY(i,j)*(EX(i,j+1)-EX(i,j) );

  for(i=1; i<N_PX-1; i++)
    for(j=1; j<N_PY-1; j++)
      Hz[ind(i,j)] = HZX(i,j) + HZY(i,j);
}
