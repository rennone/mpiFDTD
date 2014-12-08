#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "field.h"
#include "nsFdtdTM.h"
#include "models.h"
#include "function.h"

//   メンバ変数    //
static dcomplex *Ezx = NULL;
static dcomplex *Ezy = NULL;
static dcomplex *Ez = NULL; //Ez(x    ,y    ) => Ez[i,j]
static dcomplex *Hx = NULL; //Hx(x    ,y+0.5) => Hx[i,j]
static dcomplex *Hy = NULL; //Hy(x+0.5,y    ) => Hy[i,j]
static double *C_HXLY = NULL, *C_HYLX = NULL;
static double *C_EZX = NULL, *C_EZY = NULL, *C_EZXLX = NULL, *C_EZYLY = NULL;
static double *C_HX = NULL, *C_HY = NULL;
static double *EPS_EZ = NULL, *EPS_HX = NULL, *EPS_HY = NULL;

//  メンバ関数   //
static void update(void);
static void finish(void);
static void reset(void);
static void init(void);
static void calcE(void);
static void calcH(void);

static void freeMemories(void);

//:public
dcomplex* nsFdtdTM_getEz(){  return Ez;}
dcomplex* nsFdtdTM_getHx(){  return Hx;}
dcomplex* nsFdtdTM_getHy(){  return Hy;}
double* nsFdtdTM_getEpsX(){  return EPS_HX;}
double* nsFdtdTM_getEpsY(){  return EPS_HY;}
double* nsFdtdTM_getEpsZ(){  return EPS_EZ;}
double* nsFdtdTM_getEps(){  return EPS_EZ;}
void(*nsFdtdTM_getUpdate(void))(){  return update;}
void(*nsFdtdTM_getFinish(void))(){  return finish;}
void(*nsFdtdTM_getReset(void))() {  return reset;}
void(* nsFdtdTM_getInit(void))() {  return init;}


Solver* nsFdtdTM_getSolver()
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

static void update()
{
  calcH();
  calcE();  

  field_nsScatteredWaveNotUPML(Ezy, EPS_EZ, 0, 0, 1.0);
  
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  FOR_FOR(fInfo_s, i, j){
    int k = field_index(i,j);
    Ez[k] = Ezx[k] + Ezy[k];
  }
}

static bool InPML(double i, double j)
{
  int p = 0;
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  return (i < fInfo_s.N_PML+p || j < fInfo_s.N_PML+p
    || (i >= fInfo_s.N_X + fInfo_s.N_PML -p )
               || (j >= fInfo_s.N_Y + fInfo_s.N_PML -p) );
}

static void calcE()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  int dx = fInfo_s.DX, dy = fInfo_s.DY;
  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    Ezx[k] = C_EZX[k]* Ezx[k]
      + C_EZXLX[k]* ( Hy[k] - Hy[k-dx] );
  }
  
  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);
    Ezy[k] = C_EZY[k]* Ezy[k]
      - C_EZYLY[k]* ( Hx[k] - Hx[k-dy] );
  }
}

static void calcH()
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
  
  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);

    if( InPML(i,j+0.5) )
    {
      Hx[k] = C_HX[k]*Hx[k] - C_HXLY[k]* (
          Ez[k+dy] - Ez[k] );  // dy
    }
    else
    {      
      // γ/2 * dx^2dy
      dcomplex ns_operator = r_2*( (Ez[k+dy+dx] + Ez[k+dy-dx] - 2*Ez[k+dy])   
                                  -(Ez[k+dx]    + Ez[k-dx]    - 2*Ez[k]   ));      
      Hx[k] = C_HX[k]*Hx[k] - C_HXLY[k]* (
          Ez[k+dy] - Ez[k] // dy
        + ns_operator );
    }
  }
  
  FOR_FOR(fInfo_s, i, j)
  {
    int k = field_index(i,j);

    if( InPML(i+0.5,j) )
    {
      Hy[k] = C_HY[k]*Hy[k] + C_HYLX[k]*(
        Ez[k+dx] - Ez[k] );
    }
    else
    {
      
      // γ/2 * dxdy^2
      dcomplex ns_operator =  r_2*(  (Ez[k+dx+dy] + Ez[k+dx-dy] - 2*Ez[k+dx] )
                                    -(Ez[k+dy]    + Ez[k-dy]    - 2*Ez[k]    ));
      Hy[k] = C_HY[k]*Hy[k] + C_HYLX[k]*(
        Ez[k+dx] - Ez[k]//dx
        + ns_operator );     
    }
  }
}

static void finish()
{
  reset();
  freeMemories();
}

static void reset()
{
  FieldInfo fInfo = field_getFieldInfo();
  char buf[128];
  sprintf(buf, "ns_tm_%dnm.txt",fInfo.h_u_nm);
  field_outputElliptic(buf, Ez);
  memset(Hx , 0, sizeof(double complex)*N_CELL);
  memset(Hy , 0, sizeof(double complex)*N_CELL);
  memset(Ezx, 0, sizeof(double complex)*N_CELL);
  memset(Ezy, 0, sizeof(double complex)*N_CELL);
  memset(Ez , 0, sizeof(double complex)*N_CELL);
}

static void freeMemories()
{
  delete(Ez);
  delete(Ezx);
  delete(Ezy);
  delete(Hx);
  delete(Hy);
  
  delete(C_HX);
  delete(C_HY);
  delete(C_EZX);
  delete(C_EZY);

  delete(C_HXLY);
  delete(C_HYLX);
  delete(C_EZXLX);
  delete(C_EZYLY);

  delete(EPS_HX);
  delete(EPS_HY);
  delete(EPS_EZ);
}

static void allocateMemories()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  Hx  = newDComplex(fInfo_s.N_CELL);
  Hy  = newDComplex(fInfo_s.N_CELL);
  Ezy = newDComplex(fInfo_s.N_CELL);
  Ezx = newDComplex(fInfo_s.N_CELL);
  Ez  = newDComplex(fInfo_s.N_CELL);

  C_HX  = newDouble(fInfo_s.N_CELL);
  C_HY  = newDouble(fInfo_s.N_CELL);
  C_EZX = newDouble(fInfo_s.N_CELL);
  C_EZY = newDouble(fInfo_s.N_CELL);

  C_HXLY = newDouble(fInfo_s.N_CELL);
  C_HYLX = newDouble(fInfo_s.N_CELL);

  C_EZXLX = newDouble(fInfo_s.N_CELL);
  C_EZYLY = newDouble(fInfo_s.N_CELL);

  EPS_HY = newDouble(fInfo_s.N_CELL);
  EPS_HX = newDouble(fInfo_s.N_CELL);
  EPS_EZ = newDouble(fInfo_s.N_CELL);  
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
 //Ez,, Hx, Hyそれぞれでσx, σx*, σy, σy*が違う(場所が違うから)
  double sig_ez_x, sig_ez_y, sig_ez_xx, sig_ez_yy;
  double sig_hx_x, sig_hx_y, sig_hx_xx, sig_hx_yy;
  double sig_hy_x, sig_hy_y, sig_hy_xx, sig_hy_yy;

  double R = 1.0e-8;
  double M = 2.0;

  //NS用係数
  double w_s = field_getOmega();
  double k_s = field_getK();
  double u_ns = sin(w_s * 0.5) / sin(k_s * 0.5);
  printf("w_s = %lf,  k_s = %lf  u_ns = %lf\n", w_s, k_s, u_ns);
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  const double sig_max = -(M+1.0)*EPSILON_0_S*C_0_S/N_PML*log(R);
  for(int i=0; i<N_PX; i++){
    for(int j=0; j<N_PY; j++){
      int k = field_index(i,j);
      
      EPS_EZ[k] = models_eps(i,j, D_XY);
      EPS_HX[k] = models_eps(i,j+0.5, D_Y);
      EPS_HY[k] = models_eps(i+0.5,j, D_X);

      // PML領域には散乱体がないので eps = EPSILON_0_S としている.
      // (領域内はσ = 0なので関係ない)
      sig_ez_x  = sig_max*field_sigmaX(i,j);      //σ_x  for ez
      sig_ez_xx = MU_0_S/EPSILON_0_S * sig_ez_x;  //σ_x* for ez      
      sig_ez_y  = sig_max*field_sigmaY(i,j);      //σ_y  for ez
      sig_ez_yy = MU_0_S/EPSILON_0_S * sig_ez_y;  //σ_y* for ez

      sig_hx_x  = sig_max*field_sigmaX(i,j+0.5);  //σ_x  for hx
      sig_hx_xx = MU_0_S/EPSILON_0_S * sig_hx_x;  //σ_x* for hx
      sig_hx_y  = sig_max*field_sigmaY(i,j+0.5);  //σ_y  for hx
      sig_hx_yy = MU_0_S/EPSILON_0_S * sig_hx_y;  //σ_y* for hx

      sig_hy_x  = sig_max*field_sigmaX(i+0.5,j);  //σ_x  for hy
      sig_hy_xx = MU_0_S/EPSILON_0_S * sig_hy_x;  //σ_x* for hy
      sig_hy_y  = sig_max*field_sigmaY(i+0.5,j);  //σ_y  for hy
      sig_hy_yy = MU_0_S/EPSILON_0_S * sig_hy_y;  //σ_y* for hy

      // PML領域以外では α = α*, β = β* となるので,簡略化
      double a_hx_y = sig_hx_y / (2*EPSILON_0_S);
      double a_hy_x = sig_hy_x / (2*EPSILON_0_S);
      double a_ez_x = sig_ez_x / (2*EPSILON_0_S);
      double a_ez_y = sig_ez_y / (2*EPSILON_0_S);

      double b_hx_y = beta(a_hx_y);
      double b_hy_x = beta(a_hy_x);
      double b_ez_x = beta(a_ez_x);
      double b_ez_y = beta(a_ez_y);
      
      // EZ
      if( InPML(i,j) )
      {
        // Standard FDTD
        C_EZX[k]   =     field_pmlCoef(EPS_EZ[k], sig_ez_x);
        C_EZXLX[k] = field_pmlCoef_LXY(EPS_EZ[k], sig_ez_x);
        C_EZY[k]   =     field_pmlCoef(EPS_EZ[k], sig_ez_y);
        C_EZYLY[k] = field_pmlCoef_LXY(EPS_EZ[k], sig_ez_y);
      }
      else
      {
        //printf("%.2lf,  %.2lf,  %.2lf,   %.2lf  \n", a_ez_x, b_ez_x, a_ez_y, b_ez_y);
        double z_ez = sqrt(MU_0_S / EPS_EZ[k]);
        //波数kは媒質に依存する(角周波数は一定)
        double n_ez   = sqrt(EPS_EZ[k] / EPSILON_0_S); //屈折率
        double k_ez_s = k_s * n_ez;              
        double u = sin(w_s*0.5) / sin(k_ez_s*0.5);
        C_EZX[k]   = coef1(a_ez_x, b_ez_x);
        C_EZXLX[k] = u * z_ez;//coef2(a_ez_x, b_ez_x) * z_ez;
        C_EZY[k]   = coef1(a_ez_y, b_ez_y);
        C_EZYLY[k] = u * z_ez;//coef2(a_ez_y, b_ez_y) * z_ez;
      }

      // HX
      if( InPML(i, j+0.5) )
      {
        C_HX[k]    =     field_pmlCoef(MU_0_S, sig_hx_yy);
        C_HXLY[k]  = field_pmlCoef_LXY(MU_0_S, sig_hx_yy);             
      }
      else
      {
        // NonStandard FDTD
        //Δt = 1, μ(i,j) = μ0
        double z_hx = sqrt(MU_0_S / EPS_HX[k]);
        double n_hx   = sqrt(EPS_HX[k] / EPSILON_0_S);
        double k_hx_s = k_s * n_hx;      
        double u = sin(w_s*0.5) / sin(k_hx_s*0.5);
        C_HX[k]    = coef1(a_hx_y, b_hx_y);
        C_HXLY[k]  = u / z_hx;//coef2(a_hx_y, b_hx_y) /  z_hx;
      }

      // HY
      if( InPML(i+0.5, j) )
      {
        C_HY[k]    =     field_pmlCoef(MU_0_S, sig_hy_xx);
        C_HYLX[k]  = field_pmlCoef_LXY(MU_0_S, sig_hy_xx);
      }
      else
      {
        double z_hy   = sqrt(MU_0_S / EPS_HY[k]);
        double n_hy   = sqrt(EPS_HY[k] / EPSILON_0_S);
        double k_hy_s = k_s * n_hy;
        double u      = sin(w_s*0.5) / sin(k_hy_s*0.5);
        C_HY[k]       = coef1(a_hy_x, b_hy_x);
        C_HYLX[k]     = u / z_hy;//coef2(a_hy_x, b_hy_x) / z_hy;
      }

      if ( InPML(i,j) && i >= fInfo_s.N_PML && i < fInfo_s.N_PML + fInfo_s.N_X)
      {
        // printf("%.3lf   ", C_EZYLY[k]/*, i, j*/);
      }
    }

    if( i>=fInfo_s.N_PML || i<fInfo_s.N_PML + fInfo_s.N_X)
    {
      //printf("\n");
    }
  } 
}

static void init()
{
  allocateMemories();
  setCoefficient();
/*
  //計算用定数の設定
  double kx_s = k_s *          pow(2.0, -0.25);	//(2の4乗根)分の1
  double ky_s = k_s * sqrt(1 - pow(2.0, -0.5 ));
  double sin2_kx = pow(sin(kx_s/2.0), 2);
  double sin2_ky = pow(sin(ky_s/2.0), 2);
  double sin2_k  = pow(sin(k_s /2.0), 2);
  R_M =  (sin2_kx + sin2_ky - sin2_k)/(4*sin2_kx*sin2_ky); //(1-γ0)/2 = R_M, ( R_P = 1とするとおかしくなる todo?)
  R_P = 1 - R_M;*/
}
