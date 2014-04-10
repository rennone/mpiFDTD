#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#include "field.h"
#include "models.h"
#include "function.h"

int N_X;
int N_Y;
int N_CELL;
int N_PML;
int N_PX;
int N_PY;

//MPI分割したときのフィールドパラメータ
static int SUB_N_X;
static int SUB_N_Y;
static int SUB_N_PX;
static int SUB_N_PY;
static int SUB_N_CELL;
static int RANK, NPROC;
static int LTRANK, RTRANK, TPRANK,BMRANK;
static int OFFSET_X, OFFSET_Y;
static MPI_Datatype MPI_DOUBLE_COMPLEX_COLUM;

//_u : 物理量変換単位, _s:シミュレーション単位
static const int H_s = 1;
static double H_u;         //1セルの大きさ(nm)
static double time;     //ステップ数
static double ray_coef; //波をゆっくり入れる為の係数;
static double waveAngle;//入射角
static double lambda_s; //波長 
static double k_s;      //波数 
static double w_s;      //角周波数
static double T_s;      //周期 
static void (*defaultWave)(double complex* p, double* eps);
static double maxTime;

static NTFFInfo ntff_info;

static void field_setNTFF(int);
static void initMpi(void);

//:public------------------------------------//
//inline 意味ない->ヘッダに書かないと行けない
double field_toCellUnit(const double phisycalUnit)
{
  return phisycalUnit/H_u;   //セル単位に変換 
}
double field_toPhisycalUnit(const double cellUnit)
{
  return cellUnit*H_u;    //物理単位(nm)に変換
}
int field_getOffsetX()
{
  return OFFSET_X;
}
int field_getOffsetY()
{
  return OFFSET_Y;
}
int field_getSubNx()
{
  return SUB_N_X;
}
int field_getSubNy()
{
  return SUB_N_Y;
}
int field_getSubNpx()
{
  return SUB_N_PX;
}
int field_getSubNpy()
{
  return SUB_N_PY;
}
int field_getSubNcell()
{
  return SUB_N_CELL;
}

void setField(const int wid, const int hei, const double _h, const int pml, const double lambda, double maxstep)
{
  H_u = _h;
  N_PX = field_toCellUnit(wid);
  N_PY = field_toCellUnit(hei);
  N_PML = pml;  
  N_X = N_PX - 2*N_PML;
  N_Y = N_PY - 2*N_PML;
  N_CELL = N_PX * N_PY; //全セル数 
  time = 0;
  maxTime = maxstep;
  
  lambda_s = field_toCellUnit(lambda);
  k_s = 2*M_PI/lambda_s;
  w_s = C_0_S*k_s;
  T_s = 2*M_PI/w_s;

  ray_coef = 0;  
  waveAngle = 0;
  
  /* NTFF設定 */  
  ntff_info.top = N_PY - N_PML - 5;
  ntff_info.bottom = N_PML + 5;
  ntff_info.left = N_PML + 5;
  ntff_info.right = N_PX - N_PML - 5;


  double len = (ntff_info.top - ntff_info.bottom)/2;
  ntff_info.RFperC = len*2;
  ntff_info.arraySize = maxTime + 2*ntff_info.RFperC;
}

//-------------------getter-------------------//
double field_getT()
{
  return T_s;
}

double  field_getK()
{
  return k_s;
}

double  field_getRayCoef()
{
  return ray_coef;
}

double  field_getOmega()
{
  return w_s;
}

double  field_getLambda()
{
  return lambda_s;
}

double  field_getWaveAngle()
{
  return waveAngle;
}

double  field_getTime()
{
  return time;
}

double  field_getMaxTime()
{
  return maxTime;
}

NTFFInfo  field_getNTFFInfo()
{
  return ntff_info;
}

//----------------------------------------//
double field_sigmaX(const double x, const double y)
{
  const int M = 2;
  if(x<N_PML)
    return pow(1.0*(N_PML-x)/N_PML, M);
  
  else if(N_PML <= x && x < (N_X+N_PML))    
    return 0;
  
  else
    return pow(1.0*(x - (N_PX-N_PML-1))/N_PML, M);
}

double field_sigmaY(const double x, double y)
{
  const int M = 2;
  if(y<N_PML)
    return pow(1.0*(N_PML - y)/N_PML,M);
  
  else if(y>=N_PML && y<(N_Y+N_PML))
    return 0.0;

  else
    return pow(1.0*(y - (N_PY-N_PML-1))/N_PML,M);
}

//pml用の係数のひな形 Δt = 1
//ep_mu εかμ(Eの係数->ε, Hの係数-> μ
//sig  σ
double field_pmlCoef(double ep_mu, double sig)
{
  return (1.0 - sig/ep_mu)/(1.0+sig/ep_mu);
}

double field_pmlCoef_LXY(double ep_mu, double sig)
{
  return 1.0/(ep_mu + sig);
  //  return 1.0/ep_mu/(1.0 + sig/ep_mu);
}

//1次元配列に変換
 int ind(const int i, const int j)
{
  return i*N_PY + j;
}
//------------------getter-------------------------//

//------------------light method----------------------//
//点光源を返す
double complex field_pointLight(void)
{
  return ray_coef * (cos(w_s*time) + sin(w_s*time)*I);
}

//------------------light method----------------------//
 void field_nextStep(void){
  time+=1.0;
  ray_coef = 1.0 - exp(-pow(0.01*time, 2));
}

 bool field_isFinish(void){
  return time >= maxTime;
}


//---------------output method---------------//

void field_outputElliptic(const char *fileName, double complex* data)
{
  printf("output start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int ang=180; ang >=0; ang--){
    double rad = ang*M_PI/180.0;
    double x = 1.2*lambda_s*cos(rad)+N_PX/2.0;
    double y = 1.2*lambda_s*sin(rad)+N_PY/2.0;
    double norm = cnorm(cbilinear(data,x,y,N_PX,N_PY));
    fprintf(fp, "%d %lf \n", 180-ang, norm);    
  }
  
  fclose(fp);
  printf("output to %s end\n", fileName);
}

void field_outputAllDataComplex(const char *fileName, double complex *data)
{
   printf("output all data start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int i=0; i<N_PX; i++){
    for(int j=0; j<N_PY; j++){
      fprintf(fp, "%lf \n",cnorm(data[ind(i,j)]));
    }
  }
  
  fclose(fp);
  printf("output all data to %s end\n", fileName);
}


void field_outputAllDataDouble(const char *fileName, double *data)
{
   printf("output all data double start\n");
  //file open
  FILE *fp;

  if( (fp=fopen(fileName, "w") ) == NULL){
    printf("cannot open file %s \n", fileName);
    exit(2);
  }

  for(int i=0; i<N_PX; i++){
    for(int j=0; j<N_PY; j++){
      if(data[ind(i,j)] != 1.0)
        fprintf(fp, "%lf \n",data[ind(i,j)]);
    }
  }
  
  fclose(fp);
  printf("output all data to %s end\n", fileName);
}
 
static void init_mpi(void)
{
  MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
  int dim = 2;          //number of dimension is 2
  int procs[2] = {0,0};         //[0]: x方向の分割数, [1]:y方向の分割数
  int period[2] = {0,0};//境界条件, 固定境界
  MPI_Comm grid_comm;
  int reorder = 1;   //re-distribute rank flag

  MPI_Dims_create(NPROC, dim, procs);
  MPI_Cart_create(MPI_COMM_WORLD, 2, procs, period, reorder, &grid_comm);
  MPI_Cart_shift(grid_comm, 0, 1, &LTRANK, &RTRANK);
  MPI_Cart_shift(grid_comm, 1, 1, &BMRANK, &TPRANK);

  //プロセス座標において自分がどの位置に居るのか求める(何行何列に居るか)
  int coordinates[2];
  MPI_Comm_rank(grid_comm, &RANK);
  MPI_Cart_coords(grid_comm, RANK, 2, coordinates);
  
  SUB_N_X = N_PX / procs[0];
  SUB_N_Y = N_PY / procs[1];
  SUB_N_PX = SUB_N_X + 2; //のりしろの分2大きい
  SUB_N_PY = SUB_N_Y + 2; //のりしろの分2大きい
  SUB_N_CELL = SUB_N_PX*SUB_N_PY;
  
  OFFSET_X = coordinates[0] * SUB_N_X; //ランクのインデックスではなく, セル単位のオフセットなのでSUB_N_Xずれる
  OFFSET_Y = coordinates[1] * SUB_N_Y;

/*これだと, 1個のデータをSUB_N_PY跳び(次のデータまでSUB_N_PY-1個隙間がある),SUB_N_X行ぶん取ってくる事になる */
  MPI_Type_vector(SUB_N_X, 1, SUB_N_PY, MPI_C_DOUBLE_COMPLEX, &MPI_DOUBLE_COMPLEX_COLUM); 
  MPI_Type_commit(&MPI_DOUBLE_COMPLEX_COLUM);
}
