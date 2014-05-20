#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simulator.h"
#include "field.h"

#include "mpiTM_UPML.h"
#include "mpiTE_UPML.h"
#include "fdtdTM.h"
#include "fdtdTE.h"
#include "fdtdTM_upml.h"
#include "fdtdTE_upml.h"

#include "models.h"
#include <sys/time.h>

static void (*update)() = NULL;
static double complex* (*getDataX)() = NULL;
static double complex* (*getDataY)() = NULL;
static double complex* (*getDataZ)() = NULL;
static void (* finishMethod)() = NULL;
static void (* initMethod)() = NULL;
static double complex* (*getDrawData)() = NULL;
static double* (* getEpsMethod )() = NULL;

static struct timeval timer1, timer2;

static void setTM()
{
  update       = fdtdTM_getUpdate();
  initMethod   = fdtdTM_getInit();
  finishMethod = fdtdTM_getFinish();

  getEpsMethod = fdtdTM_getEps;
  getDataX = fdtdTM_getHx;
  getDataY = fdtdTM_getHy;
  getDataZ = fdtdTM_getEz;

  getDrawData = getDataZ;
  printf("TM mode \n");

  makeDirectory("TM");
  moveDirectory("TM");
}

static void setTE()
{
  update       = fdtdTE_getUpdate();
  initMethod   = fdtdTE_getInit();
  finishMethod = fdtdTE_getFinish();

  getEpsMethod = fdtdTE_getEps;
  getDataX = fdtdTE_getEx;
  getDataY = fdtdTE_getEy;
  getDataZ = fdtdTE_getHz;

  getDrawData = getDataY;
  printf("TE mode \n");
  makeDirectory("TE");
  moveDirectory("TE");
}

static void setTMupml()
{
  update       = fdtdTM_upml_getUpdate();
  initMethod   = fdtdTM_upml_getInit();
  finishMethod = fdtdTM_upml_getFinish();

  getEpsMethod = fdtdTM_upml_getEps;
  getDataX = fdtdTM_upml_getHx;
  getDataY = fdtdTM_upml_getHy;
  getDataZ = fdtdTM_upml_getEz;

  getDrawData = getDataZ;
  printf("TM UPML mode \n");
  makeDirectory("TM_UPML");
  moveDirectory("TM_UPML");
}

static void setTEupml()
{
  update       = fdtdTE_upml_getUpdate();
  initMethod   = fdtdTE_upml_getInit();
  finishMethod = fdtdTE_upml_getFinish();

  getEpsMethod = fdtdTE_upml_getEps;
  getDataX = fdtdTE_upml_getEx;
  getDataY = fdtdTE_upml_getEy;
  getDataZ = fdtdTE_upml_getHz;

  getDrawData = getDataY;
  printf("TE UPML mode \n");
  makeDirectory("TE_UPML");
  moveDirectory("TE_UPML");
}

static void setMPITMupml(){
  update       = mpi_fdtdTM_upml_getUpdate();
  initMethod   = mpi_fdtdTM_upml_getInit();
  finishMethod = mpi_fdtdTM_upml_getFinish();

  getEpsMethod = mpi_fdtdTM_upml_getEps;
  
  getDataX = mpi_fdtdTM_upml_getHx;
  getDataY = mpi_fdtdTM_upml_getHy;
  getDataZ = mpi_fdtdTM_upml_getEz;
  
  getDrawData = getDataZ;
  printf("MPI TM UPML mode \n");
  makeDirectory("MPI_TM_UPML");
  moveDirectory("MPI_TM_UPML");
}

static void setMPITEupml(){
  update       = mpi_fdtdTE_upml_getUpdate();
  initMethod   = mpi_fdtdTE_upml_getInit();
  finishMethod = mpi_fdtdTE_upml_getFinish();
  
  getEpsMethod = mpi_fdtdTE_upml_getEps;

  getDataX = mpi_fdtdTE_upml_getEx;
  getDataY = mpi_fdtdTE_upml_getEy;
  getDataZ = mpi_fdtdTE_upml_getHz;

  getDrawData = getDataY;
  
  printf("MPI TE UPML mode \n");
  makeDirectory("MPI_TE_UPML");
  moveDirectory("MPI_TE_UPML");
}

static void setSolver(enum SOLVER solver)
{
  switch(solver){
  case TM_2D:
    setTM();
    break;
  case TE_2D:
    setTE();
    break;
  case TE_UPML_2D:
    setTEupml();
    break;
  case TM_UPML_2D:
    setTMupml();
    break;
  case MPI_TM_UPML_2D:
    setMPITMupml();
    break;
  case MPI_TE_UPML_2D:
    setMPITEupml();
    break;
  default:
    printf("error, not implement simulator (simulator.c)\n");
    exit(2);
    break;
  }
  (*initMethod)(); //Solverの初期化, EPS, Coeffの設定  
}

void simulator_calc(){
  (*update)();
  
  field_nextStep();   //時間を一つ進める  
}

void simulator_init(FieldInfo field_info, enum MODEL model, enum SOLVER solver)
{
  //横幅(nm), 縦幅(nm), 1セルのサイズ(nm), pmlレイヤの数, 波長(nm), 計算ステップ
  initField(field_info);

  //NO_MODEL. MIE_CYLINDER, SHELF, NONSHELF
  // UNDONE : SHELF, NONSHELFモデル
  setModel(model);  //次にこれ,  モデル(散乱体)を定義

  setSolver(solver);//Solverの設定と初期化

  gettimeofday(&timer1, NULL); //開始時間の取得
}

void simulator_finish()
{  
  printf("finish at %d step \n", (int)field_getTime());
  gettimeofday(&timer2,NULL);
  printf("time = %lf \n", timer2.tv_sec-timer1.tv_sec+(timer2.tv_usec-timer1.tv_usec)*1e-6);

  (*finishMethod)();   //メモリの解放等, solverの終了処理
}

double complex* simulator_getDrawingData(void){
  return (* getDrawData)();
}

bool simulator_isFinish(void)
{
  return field_isFinish();
}

//屈折率のマップを取ってくる
double *simulator_getEps()
{
  return getEpsMethod();  
}
