#define USE_MATH_DEFINES
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include <sys/stat.h>
#include "mpiTM_UPML.h"
#include "field.h"
#include "models.h"
#include "function.h"
#include "ntffTM.h"

/* about MPI  */
static int rank;      //MPIのランク
static int nproc;     //全プロセス数
static int offsetX, offsetY;  //計算領域における, このミニフィールドのオフセット量
static int ltRank, rtRank, tpRank, bmRank; //左右上下のランク
static int SUB_N_X, SUB_N_Y;    //分割領域の大きさ
static int SUB_N_PX, SUB_N_PY; //のりしろ(2セル)を含めた, 分割領域の大きさ
static int SUB_N_CELL;
static MPI_Datatype X_DIRECTION_DOUBLE_COMPLEX;

/* about NTFF */
static double complex *Ux = NULL;
static double complex *Uy = NULL;
static double complex *Wz = NULL;

#ifdef DEBUG
static double complex *debug_U[4];
static double complex *debug_W[4];
#endif

static void initDebug(void);

//Ez(i    , j    ) -> Ez[i,j];
//Hx(i    , j+0.5) -> Hx[i,j];
//Hy(i+0.5, j    ) -> Hy[i,j];
static double complex *Ez = NULL;
static double complex *Jz = NULL;
static double complex *Dz = NULL;

static double complex *Hx = NULL;
static double complex *Mx = NULL;
static double complex *Bx = NULL;

static double complex *Hy = NULL;
static double complex *My = NULL;
static double complex *By = NULL;

static double *C_JZ = NULL, *C_MX = NULL, *C_MY = NULL;
static double *C_JZHXHY = NULL, *C_MXEZ = NULL, *C_MYEZ = NULL;
static double *C_DZ = NULL, *C_BX = NULL, *C_BY = NULL;
static double *C_DZJZ0 = NULL, *C_DZJZ1 = NULL;
static double *C_BXMX0 = NULL, *C_BXMX1 = NULL;
static double *C_BYMY0 = NULL, *C_BYMY1 = NULL;

static double *EPS_EZ = NULL, *EPS_HX = NULL, *EPS_HY = NULL;

static void update(void);
static void finish(void);
static void output(void);
static void freeMemories(void);
static void allocateMemories(void);

static void setCoefficient(void);
static void init(void);
static void init_mpi(void);

static inline void planeWave(double complex *p, double *eps);
static inline void scatteredWave(double complex *p, double *eps);
static inline void Connection_ISend_IRecvE(void);
static inline void Connection_ISend_IRecvH(void);
static inline void Connection_SendRecvE(void);
static inline void Connection_SendRecvH(void);
static inline void calcJD(void);
static inline void calcE(void);
static inline void calcMB(void);
static inline void calcH(void);

static inline void _CalcJD();
static inline void _CalcMB();
static inline void ntffCoef(double time, double timeShift, int *m, double *a, double *b, double *ab);
static inline void ntff(void);
static inline void ntffOutput();

//:public
inline int mpi_fdtdTM_upml_getSubNx(void)
{
  return SUB_N_X;
}
inline int mpi_fdtdTM_upml_getSubNy(void)
{
  return SUB_N_Y;
}
inline int mpi_fdtdTM_upml_getSubNpx(void)
{
  return SUB_N_PX;
}
inline int mpi_fdtdTM_upml_getSubNpy(void)
{
  return SUB_N_PY;
}
inline int mpi_fdtdTM_upml_getSubNcell(void)
{
  return SUB_N_CELL;
}

void (* mpi_fdtdTM_upml_getUpdate(void))(void)
{
  return update;
}
void (* mpi_fdtdTM_upml_getFinish(void))(void)
{
  return finish;
}
void (* mpi_fdtdTM_upml_getInit(void))(void)
{
  return init;
}
double complex* mpi_fdtdTM_upml_getHx(void){
  return Hx;
}
double complex* mpi_fdtdTM_upml_getHy(void){
  return Hy;
}
double complex* mpi_fdtdTM_upml_getEz(void){
  return Ez;
}
double* mpi_fdtdTM_upml_getEps(void){
  return EPS_EZ;
}

//:private
//--------------------getter--------------------//
static inline int subInd(const int *i, const int *j)
{
  return (*i)*SUB_N_PY + (*j);
}

static inline int subIndv(const int i, const int j)
{
  return (i)*SUB_N_PY + (j);
}

static inline int subIndLeft(const int i)
{
  return i - SUB_N_PY;
}

static inline int subIndRight(const int i)
{
  return i + SUB_N_PY;
}

static inline int subIndTop(const int i)
{
  return i+1;
}

static inline int subIndBottom(const int i)
{
  return i-1;
}

//fileOpen用
static FILE* fileOpen(const char* fileName)
{
  char name[256];
  //フォルダー以下にする.
  sprintf(name, "MPI_TM_UPML/%s", fileName);
  
  FILE *fp;
  if( (fp=fopen(name, "w") ) == NULL )
  {
    printf("cannot open file %s \n", name);
    exit(2);
  }
  return fp;
}

//Initalize
static void init(void)
{
  init_mpi();
  allocateMemories();
  setCoefficient();
  initDebug();
}

//Update
static void update(void)
{
//  _CalcJD();
  calcJD();
  calcE();

//  MPI_Barrier(MPI_COMM_WORLD); //いらない?
  scatteredWave(Ez, EPS_EZ);
  //planeWave(Ez, EPS_EZ);

  
  //Connection_SendRecvE();
  Connection_ISend_IRecvE();
  
//  _CalcMB();
  calcMB();
  calcH();
  Connection_ISend_IRecvH();
  //Connection_SendRecvH();

  ntff();
}

//Finish
static void finish(void)
{
  //最後は同期をとっておく
  MPI_Barrier(MPI_COMM_WORLD);
  
//  ntffOutput();
//  output();

  dcomplex *res = (dcomplex*)malloc(sizeof(dcomplex)*360);
  ntffTM_FrequencySplit(Hx, Hy, Ez, res);
  free(res);
  
  freeMemories();
  
  MPI_Finalize();
}

static inline void Connection_ISend_IRecvH(void)
{
  MPI_Status status;  
  MPI_Request req1, req2, req3, req4;
  //left & right
  //this needs only Hy[i-1, j] so send to right and recieve from left 
  int ltRecv = subIndv(0,1);
  int rtSend = subIndv(SUB_N_PX-2, 1);  
  MPI_Isend(&Hy[rtSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &req1);
  MPI_Irecv(&Hy[ltRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1, MPI_COMM_WORLD, &req2);
  
  //this needs only Hy[i,j-1] so send to top and recieve from bottom   
  int bmRecv = subIndv(1,0);
  int tpSend = subIndv(1,SUB_N_PY-2);
  MPI_Isend(&Hx[tpSend], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &req3);
  MPI_Irecv(&Hx[bmRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &req4);

  MPI_Wait(&req1, &status);
  MPI_Wait(&req2, &status);
  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status); 
}

static inline void Connection_ISend_IRecvE(void)
{
  MPI_Status status;
  MPI_Request req1, req2, req3, req4;

  //left & right
  //this needs only Ez[i+1, j] so send to left and recieve from right 
  int rtRecv = subIndv(SUB_N_PX-1,1);
  int ltSend = subIndv(         1,1);
  MPI_Isend(&Ez[ltSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1, MPI_COMM_WORLD, &req1);
  MPI_Irecv(&Ez[rtRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &req2);  
  //bottom & top
//this needs only Ez[i, j+1] so send to bottom and recieve from top
  int bmSend = subIndv(1,1);
  int tpRecv = subIndv(1,SUB_N_PY-1);
  MPI_Isend(&Ez[bmSend], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &req3);
  MPI_Irecv(&Ez[tpRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &req4);  
  MPI_Wait(&req1, &status);
  MPI_Wait(&req2, &status);
  MPI_Wait(&req3, &status);
  MPI_Wait(&req4, &status);
}

static inline void Connection_SendRecvE(void)
{
  MPI_Status status;
  
  //left & right
  //this needs only Ez[i+1, j] so send to left and recieve from right 
  int rtRecv = subIndv(SUB_N_PX-1,1);
  int ltSend = subIndv(         1,1);
  MPI_Sendrecv(&Ez[ltSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1,
               &Ez[rtRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1, MPI_COMM_WORLD, &status);
  
  //bottom & top
//this needs only Ez[i, j+1] so send to bottom and recieve from top
  int bmSend = subIndv(1,1);
  int tpRecv = subIndv(1,SUB_N_PY-1);
  MPI_Sendrecv(&Ez[bmSend], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1,
               &Ez[tpRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1, MPI_COMM_WORLD, &status);
}

static inline void Connection_SendRecvH(void)
{
  MPI_Status status;

  //left & right
  //this needs only Hy[i-1, j] so send to right and recieve from left 
  int ltRecv = subIndv(0,1);
  int rtSend = subIndv(SUB_N_PX-2, 1);
  MPI_Sendrecv(&Hy[rtSend], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, rtRank, 1,
               &Hy[ltRecv], SUB_N_Y, MPI_C_DOUBLE_COMPLEX, ltRank, 1,MPI_COMM_WORLD, &status);
  
  //this needs only Hy[i,j-1] so send to top and recieve from bottom   
  int bmRecv = subIndv(1,0);
  int tpSend = subIndv(1,SUB_N_PY-2);
  MPI_Sendrecv(&Hx[tpSend], 1, X_DIRECTION_DOUBLE_COMPLEX, tpRank, 1,
               &Hx[bmRecv], 1, X_DIRECTION_DOUBLE_COMPLEX, bmRank, 1, MPI_COMM_WORLD, &status);
 
}

//Standard Scattered Wave
static inline void scatteredWave(double complex *p, double *eps)
{
  double time = field_getTime();
  double w_s  = field_getOmega();
  double ray_coef = field_getRayCoef();
  double k_s = field_getK();
  double rad = field_getWaveAngle()*M_PI/180;	//ラジアン変換

//毎回計算すると時間かかりそうだから代入しておく  
  double _cos = cos(rad), _sin = sin(rad);
  double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;
  const double t0 = 100;      //t0 stepでmaxになるように時間移動
  
  //ガウシアンパルス
  const double beam_width = 50;
  
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      int k = subInd(&i,&j);
      int x = i-1+offsetX;
      int y = j-1+offsetY;

      /*
      //ガウシアンパルス
      const double r = (x*_cos+y*_sin)/C_0_S-(time-t0);
      const double gaussian_coef = exp( -pow(r/beam_width, 2 ) );
      p[k] += gaussian_coef*(EPSILON_0_S/eps[k] - 1)*cexp(I*r*w_s);
      */
      
      //単一波長の散乱波
      double kr = x*ks_cos+y*ks_sin;
      p[k] += ray_coef*(EPSILON_0_S/eps[k] - 1.0)*cexp( I*(kr-w_s*time) );
      
    }
  }
}

//平面波
static inline void planeWave(double complex *p, double *eps)
{
  const double time = field_getTime();
  const double w_s  = field_getOmega();
  const double ray_coef = field_getRayCoef();
  const double k_s = field_getK();  
  const double rad = field_getWaveAngle()*M_PI/180;	//ラジアン変換
  const double ks_cos = cos(rad)*k_s, ks_sin = sin(rad)*k_s;	//毎回計算すると時間かかりそうだから,代入しておく

  const double t0 = 500; //t0 stepでmaxになるように時間移動
  const double a = 1.0e-4;
  double gaussianCoef = exp( -M_PI*a*pow(t0 - time, 2 ) );

  NTFFInfo nInfo = field_getNTFFInfo();
  const int i=nInfo.left;
  int j;  
  for(j=1; j<SUB_N_PY-1; j++)
  {
    int k = subInd(&i,&j);
    int x = i-1+offsetX;
    int y = j-1+offsetY;
    //double kr = x*ks_cos + y*ks_sin; //k_s*(i*cos + j*sin)
    // p[k] += 100*gaussianCoef*ray_coef*cexp(I*(kr-w_s*(time-t0)));  
    double kr = (x*ks_cos + y*ks_sin) - time; //k_s*(i*cos + j*sin)
    p[k] += ray_coef*cexp(I*kr*w_s);
  }
}

//Eは別に計算した方が速い
//=>配列のアクセスの問題がfor文のインクリメントとかにかかるコストを上回っているから?
static inline void _CalcJD()
{
  int sub_py = field_getSubNpy();
  int k    = sub_py + 1;                      //i行目の0列目
  int last = field_getSubNcell() - sub_py;//最後のセル(右上)の左下のセルまで

  // for(i=1..SUB_N_PX-1)
  //   for(j=1..SUB_N_PY-1) と同じ
  for(; k<last; k+=2) {
    // 列の最後では, 次の列の1番目まで移動([i][SUB_N_PY-1] => [i+1][1])
      int endRow = k+sub_py-2;
    for( ; k<endRow; k++)
    {
      const int k_lft = subIndLeft(k);    //一つ左
      const int k_btm = subIndBottom(k);  //一つ下
      const double complex nowJz = Jz[k];
      Jz[k] = C_JZ[k]*Jz[k] + C_JZHXHY[k]*(+Hy[k] - Hy[k_lft] - Hx[k] + Hx[k_btm]);
      Dz[k] = C_DZ[k]*Dz[k] + C_DZJZ1[k]*Jz[k] - C_DZJZ0[k]*nowJz;
//      Ez[k] = Dz[k]/EPS_EZ[k];
    }
  }
}

//calculate J and D
static inline void calcJD()
{
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      const int k = subInd(&i,&j);
      const int k_lft = subIndLeft(k); //一つ左
      const int k_btm = subIndBottom(k);       //一つ下
      const double complex nowJz = Jz[k];
      Jz[k] = C_JZ[k]*Jz[k] + C_JZHXHY[k]*(+Hy[k] - Hy[k_lft] - Hx[k] + Hx[k_btm]);
      Dz[k] = C_DZ[k]*Dz[k] + C_DZJZ1[k]*Jz[k] - C_DZJZ0[k]*nowJz;
    }
  }
}

//calculate E 
static inline void calcE()
{
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      const int k = subInd(&i,&j);
      Ez[k] = Dz[k]/EPS_EZ[k];
    }
  }
}


//calculate M and B
//Hは別に計算した方が速い
static inline void _CalcMB()
{
  int sub_py = field_getSubNpy();
  int k    = sub_py + 1;                      //i行目の0列目
  int last = field_getSubNcell() - sub_py;//最後のセル(右上)の左下のセルまで

  for(; k<last; k+=2)
  {
    int end = k+sub_py-2;
    for(; k<end; k++)
    {
      const int k_top = subIndTop(k); //一つ上
      double complex nowMx = Mx[k];
      Mx[k] = C_MX[k]*Mx[k] - C_MXEZ[k]*(Ez[k_top] - Ez[k]);
      Bx[k] = C_BX[k]*Bx[k] + C_BXMX1[k]*Mx[k] - C_BXMX0[k]*nowMx;

      const int k_rht = subIndRight(k); //一つ右     
      double complex nowMy = My[k];
      My[k] = C_MY[k]*My[k] - C_MYEZ[k]*(-Ez[k_rht] + Ez[k]);
      By[k] = C_BY[k]*By[k] + C_BYMY1[k]*My[k] - C_BYMY0[k]*nowMy;

//      Hx[k] = Bx[k]/MU_0_S;
//      Hy[k] = By[k]/MU_0_S;
    }
  }  
}

//calculate M and B
static inline void calcMB()
{  
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      const int k = subInd(&i,&j);
      const int k_top = subIndTop(k); //一つ上
      double complex nowMx = Mx[k];
      Mx[k] = C_MX[k]*Mx[k] - C_MXEZ[k]*(Ez[k_top] - Ez[k]);
      Bx[k] = C_BX[k]*Bx[k] + C_BXMX1[k]*Mx[k] - C_BXMX0[k]*nowMx;
    }
  }
  
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      const int k = subInd(&i,&j);
      const int k_rht = subIndRight(k); //一つ右     
      double complex nowMy = My[k];
      My[k] = C_MY[k]*My[k] - C_MYEZ[k]*(-Ez[k_rht] + Ez[k]);
      By[k] = C_BY[k]*By[k] + C_BYMY1[k]*My[k] - C_BYMY0[k]*nowMy;
    }
  }
}

//calculate H
static inline void calcH()
{
  for(int i=1; i<SUB_N_PX-1; i++)
  {
    for(int j=1; j<SUB_N_PY-1; j++)
    {
      const int k = subInd(&i,&j);
      Hx[k] = Bx[k]/MU_0_S;
      Hy[k] = By[k]/MU_0_S;
    }
  }  
}

//----------------------------------------
//ff initialize and finaliz
//----------------------------------------
static void output()
{
  if(rank == 0){
    MPI_Status status;
    double complex *entire = (double complex*)malloc(sizeof(double complex)*N_CELL);
    memset(entire, 0, sizeof(double complex)*N_CELL);
    
    for(int i=1, x=offsetX; i<SUB_N_PX-1; i++, x++)
      for(int j=1, y=offsetY; j<SUB_N_PY-1; j++, y++)
        entire[ind(x,y)] = Ez[subInd(&i, &j)];

    int offsets[2];
    for(int i=1; i<nproc; i++)
    {
      MPI_Recv(offsets, 2, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
      MPI_Recv(Ez, SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, i, 1, MPI_COMM_WORLD, &status);
      
      for(int i=1, x=offsets[0]; i<SUB_N_PX-1; i++, x++)
        for(int j=1, y=offsets[1]; j<SUB_N_PY-1; j++, y++)
          entire[ind(x,y)] = Ez[subInd(&i, &j)];      
    }
    field_outputElliptic("mpi_mie.txt", entire);
    free(entire);
  }
  else{
    int offsets[2];
    offsets[0] = offsetX; offsets[1] = offsetY;
    MPI_Send(offsets, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
    MPI_Send(Ez, SUB_N_CELL, MPI_C_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

//free memories and MPI_TYPE
static void freeMemories(void)
{
  MPI_Type_free(&X_DIRECTION_DOUBLE_COMPLEX);
  
  if(Ez != NULL){   free(Ez); Ez = NULL;  }  
  if(Hx != NULL){   free(Hx); Hx = NULL;  }
  if(Hy != NULL){   free(Hy); Hy = NULL;  }
    
  if(C_JZ != NULL){ free(C_JZ); C_JZ = NULL;  }
  if(C_JZHXHY != NULL){ free(C_JZHXHY); C_JZHXHY = NULL;  }
  
  if(C_DZ != NULL){ free(C_DZ); C_DZ = NULL;  }
  if(C_DZJZ1 != NULL){ free(C_DZJZ1); C_DZJZ1 = NULL;  }
  if(C_DZJZ0 != NULL){ free(C_DZJZ0); C_DZJZ0 = NULL;  }

  if(C_MX != NULL){ free(C_MX); C_MX = NULL;  }
  if(C_MXEZ != NULL){ free(C_MXEZ); C_MXEZ = NULL;  }
  
  if(C_BX != NULL){ free(C_BX); C_BX = NULL;  }
  if(C_BXMX1 != NULL){ free(C_BXMX1); C_BXMX1 = NULL;  }
  if(C_BXMX0 != NULL){ free(C_BXMX0); C_BXMX0 = NULL;  }
  
  if(C_MY != NULL){ free(C_MY); C_MY = NULL;  }
  if(C_MYEZ != NULL){ free(C_MYEZ); C_MYEZ = NULL;  }

  if(C_BY != NULL){ free(C_BY); C_BY = NULL;  }
  if(C_BYMY1 != NULL){ free(C_BYMY1); C_BYMY1 = NULL;  }
  if(C_BYMY0 != NULL){ free(C_BYMY0); C_BYMY0 = NULL;  }

  if(Ux != NULL){free(Ux); Ux = NULL;}
  if(Uy != NULL){free(Uy); Uy = NULL;}
  if(Wz != NULL){free(Wz); Wz = NULL;}
}

//malloc memories and initialize with 0
static void allocateMemories(void)
{
    
  Ez = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  Dz = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  Jz = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  
  Hx = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  Mx = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  Bx = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  
  Hy = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  My = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  By = (double complex*)malloc(sizeof(double complex)*SUB_N_CELL);
  
  C_JZ = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MY = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_DZ = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BY = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_JZHXHY = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MXEZ = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_MYEZ = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_DZJZ0 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_DZJZ1 = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_BXMX1 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BXMX0 = (double *)malloc(sizeof(double)*SUB_N_CELL);

  C_BYMY1 = (double *)malloc(sizeof(double)*SUB_N_CELL);
  C_BYMY0 = (double *)malloc(sizeof(double)*SUB_N_CELL);

  EPS_HY = (double *)malloc(sizeof(double)*SUB_N_CELL);
  EPS_HX = (double *)malloc(sizeof(double)*SUB_N_CELL);
  EPS_EZ = (double *)malloc(sizeof(double)*SUB_N_CELL);

  memset(Hx, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Hy, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Ez, 0, sizeof(double complex)*SUB_N_CELL);

  memset(Mx, 0, sizeof(double complex)*SUB_N_CELL);
  memset(My, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Jz, 0, sizeof(double complex)*SUB_N_CELL);

  memset(Bx, 0, sizeof(double complex)*SUB_N_CELL);
  memset(By, 0, sizeof(double complex)*SUB_N_CELL);
  memset(Dz, 0, sizeof(double complex)*SUB_N_CELL);

  NTFFInfo info = field_getNTFFInfo();
  Ux = (double complex*)malloc(sizeof(double complex)*360*info.arraySize);
  Uy = (double complex*)malloc(sizeof(double complex)*360*info.arraySize);
  Wz = (double complex*)malloc(sizeof(double complex)*360*info.arraySize);
  memset(Ux, 0, sizeof(double complex)*360*info.arraySize);
  memset(Uy, 0, sizeof(double complex)*360*info.arraySize);
  memset(Wz, 0, sizeof(double complex)*360*info.arraySize);
}

static void setCoefficient(void)
{
  //Ez,, Hx, Hyそれぞれでσx,σyが違う(場所が違うから)
  double sig_ez_x, sig_ez_y;
  double sig_hx_x, sig_hx_y;
  double sig_hy_x, sig_hy_y;  
  const double R = 1.0e-8;
  const double M = 2.0;
  //const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML*log(R);
  const double sig_max = -(M+1.0)*EPSILON_0_S*LIGHT_SPEED_S/2.0/N_PML/cos(M_PI/3)*log(R);
  
  for(int i=1; i<SUB_N_PX-1; i++){
    for(int j=1; j<SUB_N_PY-1; j++){
      int k = subInd(&i,&j);
      int x = i-1 + offsetX; //1から始めているので-1してある(0は隣の領域の値が入るから)
      int y = j-1 + offsetY;
      
      EPS_EZ[k] = models_eps(x,y, D_XY);     
      EPS_HX[k] = models_eps(x,y+0.5, D_Y);
      EPS_HY[k] = models_eps(x+0.5,y, D_X);

      sig_ez_x = sig_max*field_sigmaX(x,y);
      sig_ez_y = sig_max*field_sigmaY(x,y);

      sig_hx_x = sig_max*field_sigmaX(x,y+0.5);
      sig_hx_y = sig_max*field_sigmaY(x,y+0.5);
      
      sig_hy_x = sig_max*field_sigmaX(x+0.5,y);
      sig_hy_y = sig_max*field_sigmaY(x+0.5,y);

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

static void init_mpi(void)
{
  MPI_Comm_size(MPI_COMM_WORLD, &nproc); //プロセス数の取得
  int dim = 2;          //number of dimension
  int procs[2] = {0,0}; //[0]: x方向の分割数, [1]:y方向の分割数 がはいる
  int period[2] = {0,0};//境界条件, 0は固定境界
  MPI_Comm grid_comm;
  int reorder = 1;   //re-distribute rank flag

  MPI_Dims_create(nproc, dim, procs); //縦横を何分割にするか自動的に計算
  MPI_Cart_create(MPI_COMM_WORLD, 2, procs, period, reorder, &grid_comm); //領域を自動分割 => procs, grid_commは変更される
  MPI_Cart_shift(grid_comm, 0, 1, &ltRank, &rtRank);
  MPI_Cart_shift(grid_comm, 1, 1, &bmRank, &tpRank);

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[2];
  MPI_Comm_rank(grid_comm, &rank);
  MPI_Cart_coords(grid_comm, rank, 2, coordinates);
  
  SUB_N_X = N_PX / procs[0];
  SUB_N_Y = N_PY / procs[1];
  SUB_N_PX = SUB_N_X + 2; //のりしろ(となりの領域の値が入る部分)の分2大きい
  SUB_N_PY = SUB_N_Y + 2;
  SUB_N_CELL = SUB_N_PX*SUB_N_PY;  
  offsetX = coordinates[0] * SUB_N_X; //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  offsetY = coordinates[1] * SUB_N_Y;

/*これだと, 1個のデータをSUB_N_PY跳び(次のデータまでSUB_N_PY-1個隙間がある),SUB_N_X行ぶん取ってくる事になる */
  MPI_Type_vector(SUB_N_X, 1, SUB_N_PY, MPI_C_DOUBLE_COMPLEX, &X_DIRECTION_DOUBLE_COMPLEX); 
  MPI_Type_commit(&X_DIRECTION_DOUBLE_COMPLEX);
}

//----------------------------------------
//ff NTFF method
//----------------------------------------
//ntff用のフォーマットでdataをfileNameに保存
static inline void ntffSaveData(const char *fileName, double complex *data)
{
  char realFile[256], imagFile[256];
  FILE *fpR, *fpI;
  NTFFInfo nInfo = field_getNTFFInfo();
  const int maxTime = field_getMaxTime();
  int ang;
    
  sprintf(realFile, "%s_r.txt", fileName);
  sprintf(imagFile, "%s_i.txt", fileName);
  
  fpR = fileOpen(realFile);
  fpI = fileOpen(imagFile);
  
  for(ang=0; ang<360; ang++)
  {
    int k= ang*nInfo.arraySize;
    int i;
    for(i=0; i < maxTime; i++)
    {
        fprintf(fpR,"%.20lf " , creal(data[k+i]));
        fprintf(fpI,"%.20lf " , cimag(data[k+i]));  
    }
    fprintf(fpR,"\n");
    fprintf(fpI,"\n");
  }
  printf("saved at %s & %s\n",realFile, imagFile);
}

static inline void debug_ntffOutput()
{
  #ifdef DEBUG
  char fileName[256];
  int i;
  for(i=0; i<4;i++)
  {
    sprintf(fileName, "debug%d_U",i);
    ntffSaveData(fileName, debug_U[i]);
    
    sprintf(fileName, "debug%d_W",i);
    ntffSaveData(fileName, debug_W[i]);
  }
  #endif
}

static inline void ntffOutput()
{
  const double w_s = field_getOmega();
  const double complex coef = csqrt( 2*M_PI*C_0_S/(I*w_s) );
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();
  double complex *Eth = (double complex*)malloc(sizeof(double complex)*360*nInfo.arraySize);
  double complex *Eph = (double complex*)malloc(sizeof(double complex)*360*nInfo.arraySize);
  int ang;
  double theta = 0;
  for(ang=0; ang<360; ang++)
  {
    double phi = ang*M_PI/180.0;
    
    int k= ang*nInfo.arraySize;    
    double sx = cos(theta)*cos(phi);
    double sy = cos(theta)*sin(phi);
    double sz = -cos(theta); //宇野先生の本では -sin(theta)になってる
    double px = -sin(phi);
    double py = cos(phi);
    int i;
    for(i=0; i < maxTime; i++)
    {
      double complex WTH = 0 + 0 + Wz[k+i]*sz;
      double complex WPH = 0 + 0;
      double complex UTH = Ux[k+i]*sx + Uy[k+i]*sy + 0;
      double complex UPH = Ux[k+i]*px + Uy[k+i]*py;      
      double complex ETH  = coef*(-Z_0_S*WTH-UPH);
      double complex EPH  = coef*(-Z_0_S*WPH+UTH);
      
      Eth[k+i] = ETH;
      Eph[k+i] = EPH;
    }
  }
  ntffSaveData("Eph", Eph);
  ntffSaveData("Eth", Eth);
  free(Eph);
  free(Eth);
  //debug_ntffOutput();
}

static inline void ntffCoef(double time, double timeShift, int *m, double *a, double *b, double *ab)
{
  double t = time + timeShift;
  *m = floor(t + 0.5); //四捨五入
  *a = (0.5 + t - *m);
  *b = 1.0-*a;
  *ab = *a-*b;
}

static void ntff()
{  
  NTFFInfo nInfo = field_getNTFFInfo();
  const double C = LIGHT_SPEED_S;
  const double R = 1.0e6;
  const double coef = 1.0/(4*M_PI*C*R);

  //分割領域系に置ける, 中心の位置
  const int cx = N_PX/2 - offsetX;
  const int cy = N_PY/2 - offsetY;
  
  int m_e, m_h;
  double a_e, b_e, ab_e, a_h, b_h, ab_h;  
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2

  //分割領域における積分路, 有効範囲は(1~SUB_N_PX-2)まで
  int tp =    nInfo.top - offsetY; //上面
  int bm = nInfo.bottom - offsetY; //下面
  int rt = nInfo.right  - offsetX; //右
  int lt = nInfo.left   - offsetX; //左

  int ang;
  //360°方向の, 遠方界を求める
  for(ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);

    const int stp = ang*nInfo.arraySize;  //角度ごとのインデックス
    //bottom side
    //normal vector n is (0,-1,0)
    //W = Js = n × H = ( 0, 0, Hx)  U = Ms = E × n = (Ez, 0,  0)
    if (  0 < bm && bm < SUB_N_PY-1)
    {
      //積分路の端で無ければ, 分割領域の端から端までが積分路である
      // left <= x < right という繰り返し条件のため, 1 ~ SUB_N_PXとなっている 
      const int subLeft  = max(1, lt);
      const int subRight = min(SUB_N_PX, rt);
      for(int i=subLeft; i<subRight; i++)
      {
        //原点との距離
        const double r2x = i  - cx;
        const double r2y = bm - cy;
        
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&i, &bm);
        const double complex ez = Ez[k];
        const double complex hx = 0.5 * ( Hx[k] + Hx[subIndBottom(k)] );
              
        Ux[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hx*b_h*coef;
        Ux[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hx*ab_h*coef;
        Ux[stp+m_e+1] -= ez*a_e*coef;        
        Wz[stp+m_h+1] -= hx*a_h*coef;
        
#ifdef DEBUG
        debug_U[0][stp+m_e-1] += ez*b_e*coef;
        debug_W[0][stp+m_h-1] += hx*b_h*coef;
        debug_U[0][stp+m_e]   += ez*ab_e*coef;
        debug_W[0][stp+m_h]   += hx*ab_h*coef;
        debug_U[0][stp+m_e+1] -= ez*a_e*coef;        
        debug_W[0][stp+m_h+1] -= hx*a_h*coef;
#endif
      }
    }
    
    //right side
    //normal vector n is (1,0,0)
    //Js = n × H = (0, 0,Hy)  Ms = E × n = ( 0,Ez,0)
    //Js -> W                 Ms -> U
    if ( 0 < rt && rt < SUB_N_PX-1)
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ )
      {
        const double r2x = rt-cx;
        const double r2y =  j-cy;

        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&rt, &j);
        const double complex ez = Ez[k];
        const double complex hy = 0.5 * ( Hy[k] + Hy[subIndLeft(k)] );      
        
        Uy[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hy*b_h*coef;      
        Uy[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hy*ab_h*coef;      
        Uy[stp+m_e+1] -= ez*a_e*coef;              
        Wz[stp+m_h+1] -= hy*a_h*coef;

#ifdef DEBUG
        debug_U[1][stp+m_e-1] += ez*b_e*coef;
        debug_W[1][stp+m_h-1] += hy*b_h*coef;
        debug_U[1][stp+m_e]   += ez*ab_e*coef;
        debug_W[1][stp+m_h]   += hy*ab_h*coef;
        debug_U[1][stp+m_e+1] -= ez*a_e*coef;        
        debug_W[1][stp+m_h+1] -= hy*a_h*coef;
#endif

      }
    }

        //top side
    //normal vector n is (0,1,0)
    //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    //Js -> W                  Ms -> U
    if ( 0 < tp && tp < SUB_N_PY-1)
    {      
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt);
      for ( int i=subLeft; i<subRight; i++ )
      {
        const double r2x  =  i-cx;
        const double r2y  = tp-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC; 
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&i, &tp);
        const double complex ez = -Ez[k];
        const double complex hx = -0.5 * ( Hx[k] + Hx[subIndBottom(k)] );

        Ux[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hx*b_h*coef;
        Ux[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hx*ab_h*coef;
        Ux[stp+m_e+1] -= ez*a_e*coef;        
        Wz[stp+m_h+1] -= hx*a_h*coef;

#ifdef DEBUG
        debug_U[2][stp+m_e-1] += ez*b_e*coef;
        debug_W[2][stp+m_h-1] += hx*b_h*coef;
        debug_U[2][stp+m_e]   += ez*ab_e*coef;
        debug_W[2][stp+m_h]   += hx*ab_h*coef;
        debug_U[2][stp+m_e+1] -= ez*a_e*coef;        
        debug_W[2][stp+m_h+1] -= hx*a_h*coef;
#endif
      }
    }

    //left side
    //normal vector n is (-1,0,0)
    //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    //Js -> W                   Ms -> U
    // (left,top) -> (left,bottom)
    if ( 0 < lt && lt < SUB_N_PX )
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ )
      {
        const double r2x = lt-cx;
        const double r2y =  j-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&lt, &j);
        const double complex ez = -Ez[k];
        const double complex hy = -0.5 * ( Hy[k] + Hy[subIndLeft(k)] );
      
        Uy[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hy*b_h*coef;      
        Uy[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hy*ab_h*coef;      
        Uy[stp+m_e+1] -= ez*a_e*coef;      
        Wz[stp+m_h+1] -= hy*a_h*coef;

#ifdef DEBUG
        debug_U[3][stp+m_e-1] += ez*b_e*coef;
        debug_U[3][stp+m_e]   += ez*ab_e*coef;
        debug_U[3][stp+m_e+1] -= ez*a_e*coef;
        debug_W[3][stp+m_h-1] += hy*b_h*coef;
        debug_W[3][stp+m_h]   += hy*ab_h*coef;
        debug_W[3][stp+m_h+1] -= hy*a_h*coef;
#endif
      }
    }
  }
}




//--------------------------------------------------//
// ntff.h ntff.c に移行済み
//--------------------------------------------------//
/*
static inline void ntffFrequency()
{
    //分割領域系に置ける, 中心の位置
  const double cx = N_PX/2 - offsetX;
  const double cy = N_PY/2 - offsetY;

  double k_s = field_getK();

  double r0 = 1.0e6;
  double complex coef = csqrt( I*k_s/(8*M_PI*r0) ) * cexp(I*k_s*r0);
  
  NTFFInfo nInfo = field_getNTFFInfo();
  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    nInfo.top - offsetY; //上面
  int bm = nInfo.bottom - offsetY; //下面
  int rt = nInfo.right  - offsetX; //右
  int lt = nInfo.left   - offsetX; //左

  const int max_angle = 360;
  int ang;
  static double complex ntffEz[360];

  printf("ntffFrequency start\n");
  for(ang=0; ang<max_angle; ang++)
  {
    double rad = ang*M_PI/180.0;
    double rx  = cos(rad), ry = sin(rad);
    double r2x, r2y;

    double complex Nz = 0;
    double complex Lx = 0;
    double complex Ly = 0;
    double complex C_EZ, C_HX, C_HY;
		
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    //積分路が分割領域内にあるか確認
    if (  0 < bm && bm < SUB_N_PY-1)
    {
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで
      for ( int i=subLeft; i<subRight; i++ )
      {
        r2x  =  i-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        r2y  = bm-cy;

        int k = subInd(&i, &bm);
        C_EZ = Ez[k];
        C_HX = 0.5 * ( Hx[k] + Hx[subIndBottom(k)] );
        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0   ,0,C_HX)
        Lx  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (C_EZ,0,   0)
      }
    }
    
    // (right,bottom) -> (right,top) n=(1,0)
    if ( 0 < rt && rt < SUB_N_PX-1)
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
       for ( int j=subBottom; j<subTop; j++ )
      {
        r2x  = rt-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        r2y  =  j-cy;

        int k = subInd(&rt, &j);
        C_EZ = Ez[k];
        C_HY = 0.5 * ( Hy[k] + Hy[subIndLeft(k)] );

        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,   0, C_HY)
        Ly  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0,C_EZ,    0)
      }
    }

    // (right,top) -> (left,top)  n=(0,1)
    if ( 0 < tp && tp < SUB_N_PY-1)
    {      
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで      
      for ( int i=subLeft; i<subRight; i++ )
      {
        r2x  =  i-cx;
        r2y  = tp-cy;

        int k = subInd(&i, &tp);
        C_EZ = Ez[k];
        C_HX = 0.5 * ( Hx[k] + Hx[subIndBottom(k)] );

        double innerProd = rx*r2x  + ry*r2y;  //内積
        Nz   -= C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0,    0, -C_HX)
        Lx   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
      }
    }

    // (left,top) -> (left,bottom)
    if ( 0 < lt && lt < SUB_N_PX )
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm); 
      for ( int j=subBottom; j<subTop; j++ )
      {
        r2x =   lt-cx; r2y = j-cy;

        int k = subInd(&lt, &j);
        C_EZ  = Ez[k];
        C_HY  = 0.5 * ( Hy[k] + Hy[subIndLeft(k)] );

        double innerProd = rx*r2x  + ry*r2y;  //内積      
        Nz   -= C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
        Ly   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
      }
    }

    double complex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    ntffEz[ang] = coef * ( Z_0_S*Nz + Lphi );
  }

  const char nameRe[256]="ntffRe.txt";
  const char nameIm[256]="ntffIm.txt";
  FILE *fpR = fileOpen(nameRe);
  FILE *fpI = fileOpen(nameIm);
  
  for(ang = 0; ang<max_angle; ang++)
  {
    fprintf(fpR, "%.18lf\n", creal(ntffEz[ang]));
    fprintf(fpI, "%.18lf\n", cimag(ntffEz[ang]));
  }

  fclose(fpR);
  fclose(fpI);
  printf("saved at %s & %s \n", nameRe, nameIm);
  printf("ntffFrequency end\n");
}
*/
//============================================================
// for Debug 
//============================================================
static void initDebug()
{
#ifdef DEBUG
  printf("debug mode\n");
  NTFFInfo info = field_getNTFFInfo();
  int i;
  for(i=0; i<4; i++)
  {    
    debug_U[i] = (double complex*)malloc(sizeof(double complex)*360*info.arraySize);
    debug_W[i] = (double complex*)malloc(sizeof(double complex)*360*info.arraySize);
  }
#endif  
}
