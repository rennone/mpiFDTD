#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nsFdtdTE.h"
#include "field.h"
#include "models.h"
#include "function.h"

//   メンバ変数    //
static dcomplex *Ex = NULL;
static dcomplex *Ey = NULL;
static dcomplex *Hz = NULL;
static dcomplex *Hzx = NULL;
static dcomplex *Hzy = NULL;
static double *C_EX = NULL, *C_EY = NULL, *C_EXLY = NULL, *C_EYLX = NULL;
static double *C_HZX= NULL, *C_HZY= NULL, *C_HZXLX= NULL, *C_HZYLY=NULL;
static double *EPS_EX=NULL, *EPS_EY=NULL, *EPS_HZ=NULL;

//------プロトタイプ宣言--------//
static void update(void);
static void finish(void);
static void reset(void);
static void init(void);
static inline void calcE(void);
static inline void calcH(void);


//:public
dcomplex* nsFdtdTE_getHz(){  return Hz;}
dcomplex* nsFdtdTE_getEx(){  return Ex;}
dcomplex* nsFdtdTE_getEy(){  return Ey;}
double* nsFdtdTE_getEpsX(){  return EPS_EX;}
double* nsFdtdTE_getEpsY(){  return EPS_EY;}
double* nsFdtdTE_getEpsZ(){  return EPS_HZ;}
double* nsFdtdTE_getEps(){  return EPS_EY;}
void(*nsFdtdTE_getUpdate(void))(){  return update;}
void(*nsFdtdTE_getFinish(void))(){  return finish;}
void(*nsFdtdTE_getReset(void))() {  return reset;}
void(* nsFdtdTE_getInit(void))() {  return init;}

Solver* nsFdtdTE_getSolver()
{
  static Solver solver;
  /*
  solver.update   = update;
  solver.finish   = finish;
  solver.init     = init;
  solver.reset    = reset;
  solver.getDataX = getHx;
  solver.getDataY = getHy;
  solver.getDataZ = getEz;
  solver.getEpsX  = getEpsX;
  solver.getEpsY  = getEpsY;
  solver.getEpsZ  = getEpsZ;
*/  
  return &solver;
}

//:private
#define FOR_FOR(fInfo_s, i, j) \
  for(int i=1; i<fInfo_s.N_PX-1; i++)          \
    for(int j=1; j<fInfo_s.N_PY-1; j++)

//-----------------領域の設定とメモリ確保-------------//
static void allocateMemories()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  Ex  = newDComplex(fInfo_s.N_CELL);   //Ex(i+0.5,j) -> Ex[i,j]
  Ey  = newDComplex(fInfo_s.N_CELL);   //Ey(i,j+0.5) -> Ey[i,j]
  Hzy = newDComplex(fInfo_s.N_CELL);  //Hz(i+0.5, j+0.5)->Hz[i,j]
  Hzx = newDComplex(fInfo_s.N_CELL);
  Hz  = newDComplex(fInfo_s.N_CELL);

  C_EX  = newDouble(fInfo_s.N_CELL);
  C_EY  = newDouble(fInfo_s.N_CELL);
  C_EXLY = newDouble(fInfo_s.N_CELL);
  C_EYLX = newDouble(fInfo_s.N_CELL);
  C_HZY = newDouble(fInfo_s.N_CELL);
  C_HZX = newDouble(fInfo_s.N_CELL);
  C_HZXLX = newDouble(fInfo_s.N_CELL);
  C_HZYLY = newDouble(fInfo_s.N_CELL);

  EPS_EY = newDouble(fInfo_s.N_CELL);
  EPS_EX = newDouble(fInfo_s.N_CELL);
  EPS_HZ = newDouble(fInfo_s.N_CELL);  
}

static bool InPML(double i, double j)
{
  int p = 0;
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  return (i < fInfo_s.N_PML+p || j < fInfo_s.N_PML+p
          || (i >= fInfo_s.N_X + fInfo_s.N_PML -p )
          || (j >= fInfo_s.N_Y + fInfo_s.N_PML -p) );
}

//簡略化したUns2
static double Uns(double alpha, double beta)
{
  double w_s = field_getOmega();
  double k_s = field_getK();
  double sinh_2_a = pow( sinh(alpha) , 2);
  double sin_2_w  = pow( sin(w_s*0.5), 2);
  double cosh_a   = cosh(2*alpha);


  return sqrt(( (sinh_2_a + sin_2_w)/cosh(alpha) - beta*beta ))  / sin(k_s * 0.5);
}

//PML用のα
static double alpha(double sigma, double ep_mu)
{
  return sigma / (2 *ep_mu);
}

//PML用に簡略化したβ
static double beta(double alpha)
{  
  return tanh(alpha) / (1 + pow(tanh(alpha),2) );
}

static double coef1(double alpha, double beta)
{
  return (1 - beta) / (1 + beta);
}

static double coef2(double alpha, double beta)
{
  double u_ns = Uns(alpha, beta);
  return u_ns / (1+beta);
}

static void setCoefficient()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  
  //Hz, Ex, Eyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ex_x, sig_ex_xx, sig_ex_y, sig_ex_yy;
  double sig_ey_x, sig_ey_xx, sig_ey_y, sig_ey_yy;
  double sig_hz_x, sig_hz_xx, sig_hz_y, sig_hz_yy;
  double R = 1.0e-8;
  double M = 2.0;
  const double sig_max = -(M+1.0)*EPSILON_0_S*C_0_S/fInfo_s.N_PML*log(R);  
  
  //NS用係数
  double w_s = field_getOmega();
  double k_s = field_getK();  

  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    EPS_EX[k] = models_eps(i+0.5,j, D_Y);
    EPS_EY[k] = models_eps(i,j+0.5, D_X);
    EPS_HZ[k] = 0.5*(models_eps(i+0.5,j+0.5, D_X) + models_eps(i+0.5,j+0.5, D_Y));

    // (領域内はσ = 0なので関係ない)
    sig_ex_x  = sig_max*field_sigmaX(i+0.5,j);  //σ_x for ex
    sig_ex_xx = MU_0_S/EPSILON_0_S * sig_ex_x;  //σ_x* for ex
    sig_ex_y  = sig_max*field_sigmaY(i+0.5,j);  //σ_y for ex
    sig_ex_yy = MU_0_S/EPSILON_0_S * sig_ex_y;  //σ_y* for ex

    sig_ey_x  = sig_max*field_sigmaX(i,j+0.5);  //σ_x  for ey
    sig_ey_xx = MU_0_S/EPSILON_0_S * sig_ey_x;  //σ_x* for ey
    sig_ey_y  = sig_max*field_sigmaY(i,j+0.5);  //σ_y  for ey
    sig_ey_yy = MU_0_S/EPSILON_0_S * sig_ey_y;  //σ_y* for ey

    sig_hz_x = sig_max*field_sigmaX(i+0.5,j+0.5); //σ_x  for hz
    sig_hz_xx = MU_0_S/EPSILON_0_S * sig_hz_x;    //σ_x* for hz
    sig_hz_y = sig_max*field_sigmaY(i+0.5,j+0.5); //σ_y for hz
    sig_hz_yy = MU_0_S/EPSILON_0_S * sig_hz_y;    //σ_y* for hz

    // PML領域以外では α = α*, β = β* となるので,簡略化
    double a_ex_y = sig_ex_y / (2*EPS_EX[k]);
    double a_ey_x = sig_ey_x / (2*EPS_EY[k]);
    double a_hz_x = sig_hz_x / (2*EPS_HZ[k]);
    double a_hz_y = sig_hz_y / (2*EPS_HZ[k]);

    double b_ex_y = beta(a_ex_y);
    double b_ey_x = beta(a_ey_x);
    double b_hz_x = beta(a_hz_x);
    double b_hz_y = beta(a_hz_y);
      
    // HZ
    if( InPML(i+0.5,j+0.5) )
    {
      // Standard FDTD
      C_HZX[k]   =     field_pmlCoef(MU_0_S, sig_hz_xx);
      C_HZXLX[k] = field_pmlCoef_LXY(MU_0_S, sig_hz_xx);
      C_HZY[k]   =     field_pmlCoef(MU_0_S, sig_hz_yy);
      C_HZYLY[k] = field_pmlCoef_LXY(MU_0_S, sig_hz_yy);
    }
    else
    {
      double z_hz = sqrt(MU_0_S / EPS_HZ[k]);
      //波数kは媒質に依存する(角周波数は一定)
      double n_hz   = sqrt(EPS_HZ[k] / EPSILON_0_S); //屈折率
      double k_hz_s = k_s * n_hz;              
      double u = sin(w_s*0.5) / sin(k_hz_s*0.5);
      C_HZX[k]   = coef1(a_hz_x, b_hz_x);
      C_HZXLX[k] = u / z_hz;//coef2(a_ez_x, b_ez_x) * z_ez;
      C_HZY[k]   = coef1(a_hz_y, b_hz_y);
      C_HZYLY[k] = u / z_hz;//coef2(a_ez_y, b_ez_y) * z_ez;
    }

    // EX
    if( InPML(i+0.5, j) )
    {
      C_EX[k]    =     field_pmlCoef(EPS_EX[k], sig_ex_y);
      C_EXLY[k]  = field_pmlCoef_LXY(EPS_EX[k], sig_ex_y);             
    }
    else
    {
      // NonStandard FDTD
      //Δt = 1, μ(i,j) = μ0
      double z_ex   = sqrt(MU_0_S / EPS_EX[k]);
      double n_ex   = sqrt(EPS_EX[k] / EPSILON_0_S);
      double k_ex_s = k_s * n_ex;      
      double u = sin(w_s*0.5) / sin(k_ex_s*0.5);
      C_EX[k]    = coef1(a_ex_y, b_ex_y);
      C_EXLY[k]  = u * z_ex;//coef2(a_hx_y, b_hx_y) /  z_hx;
    }

    // EY
    if( InPML(i, j+0.5) )
    {
      C_EY[k]    =     field_pmlCoef(EPS_EY[k], sig_ey_x);
      C_EYLX[k]  = field_pmlCoef_LXY(EPS_EY[k], sig_ey_x);
    }
    else
    {
      double z_ey   = sqrt(MU_0_S / EPS_EY[k]);
      double n_ey   = sqrt(EPS_EY[k] / EPSILON_0_S);
      double k_ey_s = k_s * n_ey;
      double u      = sin(w_s*0.5) / sin(k_ey_s*0.5);
      C_EY[k]       = coef1(a_ey_x, b_ey_x);
      C_EYLX[k]     = u * z_ey;//coef2(a_hy_x, b_hy_x) / z_hy;
    }    
  } 
}

static void init()
{
  allocateMemories();
  setCoefficient();  
}

static void freeMemories()
{
  delete(Hz);
  delete(Hzx);
  delete(Hzy);
  delete(Ex);
  delete(Ey);
  
  delete(C_EX);
  delete(C_EY);
  delete(C_HZX);
  delete(C_HZY);

  delete(C_EXLY);
  delete(C_EYLX);
  delete(C_HZXLX);
  delete(C_HZYLY);

  delete(EPS_EX);
  delete(EPS_EY);
  delete(EPS_HZ);
}
//---------------------メモリの解放--------------------//
static void finish()
{
  reset();
  freeMemories();
}

static void reset()
{
  FieldInfo fInfo = field_getFieldInfo();
  char buf[128];
  sprintf(buf, "ns_te_%dnm.txt",fInfo.h_u_nm);
  field_outputElliptic(buf, Ey);

  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  memset(Ex, 0, sizeof(dcomplex)*fInfo_s.N_CELL);
  memset(Ey, 0, sizeof(dcomplex)*fInfo_s.N_CELL);
  memset(Hzx,0, sizeof(dcomplex)*fInfo_s.N_CELL);
  memset(Hzy,0, sizeof(dcomplex)*fInfo_s.N_CELL);
  memset(Hz ,0, sizeof(dcomplex)*fInfo_s.N_CELL);
}

static void update(){
  calcH();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  FOR_FOR(fInfo_s, i, j){
    int k = field_index(i,j);
    Hz[k] = Hzx[k] + Hzy[k];
  }
  
  calcE();

  WaveInfo_S wInfo = field_getWaveInfo_S();
  double co = cos( (wInfo.Angle_deg+90) * M_PI/ 180.0);
  double si = sin( (wInfo.Angle_deg+90) * M_PI/ 180.0);

  if(co != 0.0)
    field_nsScatteredWaveNotUPML(Ex, EPS_EX, 0, 0.5, co);
  if(si != 0.0)
    field_nsScatteredWaveNotUPML(Ey, EPS_EY, 0, 0.5, si);
}

//電界の計算 
static inline void calcE(void)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  int dx = fInfo_s.DX, dy = fInfo_s.DY;

  double k_s = field_getK();
  double r =  1.0/6.0 + k_s*k_s/180.0 - pow(k_s,4) / 23040;
  double r_2 = r/2.0;

  double kx_s = k_s * pow(2.0, -0.25);
  double ky_s = k_s * sqrt(1 - pow(2.0, -0.5));
  double sin2_kx = pow(sin(kx_s*0.5), 2);
  double sin2_ky = pow(sin(ky_s*0.5), 2);
  double sin2_k  = pow(sin(k_s*0.5), 2);

  double Rm = (sin2_kx + sin2_ky - sin2_k) / (4*sin2_kx*sin2_ky);
  double Rp = 1 - Rm;

  //Ex
  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    if( InPML(i+0.5, j) )
    {
      Ex[k] = C_EX[k]*Ex[k] + C_EXLY[k]*( Hz[k] - Hz[k-dy]);
    }
    else
    {
      dcomplex ns_operator = r_2*(  (Hz[k+dx] + Hz[k-dx] - 2*Hz[k])
                                   -(Hz[k-dy+dx] + Hz[k-dy-dx] - 2*Hz[k-dy]));
      Ex[k] = C_EX[k]*Ex[k] + C_EXLY[k] * ( Hz[k] - Hz[k-dy]
                                            + ns_operator  );
    }
  }

  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    if( InPML(i, j+0.5) )
    {
      Ey[k] = C_EY[k]*Ey[k] - C_EYLX[k]*( Hz[k] - Hz[k-dx]);
    }
    else
    {
      dcomplex ns_operator = r_2*(  (Hz[k+dy] + Hz[k-dy] - 2*Hz[k])
                                   -(Hz[k-dx+dy] + Hz[k-dx-dy] - 2*Hz[k-dx]) );
      Ey[k] = C_EY[k]*Ey[k] - C_EYLX[k]*( Hz[k] - Hz[k-dx]
                                          + ns_operator);
    }
  }
}

//磁界の計算 
static inline void calcH()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  int dx = fInfo_s.DX, dy = fInfo_s.DY;
  FOR_FOR(fInfo_s, i, j)
  {
    int k  = field_index(i,j);
    Hzx[k] = C_HZX[k]*Hzx[k] - C_HZXLX[k]*(Ey[k+dx]-Ey[k]); 
  }

  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    Hzy[k] = C_HZY[k]*Hzy[k] + C_HZYLY[k]*(Ex[k+dy]-Ex[k]);
  }
}
