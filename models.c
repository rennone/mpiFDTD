//#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
#include "concentricCircleModel.h"
#include "multiLayerModel.h"
#include "morphoScaleModel.h"
#include "bool.h"
#include "function.h"

static void (*initModelMethod)();
static void (*needSizeMethod)(int*, int*);
static double (*epsMethod)(double, double, int, int);
static bool (*isFinishMethod)(void);
static void (*moveDirectoryMethod)(void);
static char *dir;

static void noModel(void)
{
  //no material
  dir = "NoModel";
  epsMethod = noModel_EPS();
  isFinishMethod = noModel_isFinish;
  printf("not implemented no Model");
  exit(2);
}

static void circleModel(void)
{
  dir = "MieCylinderModel";
  epsMethod           = circleModel_EPS();
  isFinishMethod      = circleModel_isFinish;  
  needSizeMethod      = circleModel_needSize;
  initModelMethod     = circleModel_init;
  moveDirectoryMethod = circleModel_moveDirectory;  
}

static void multiLayerModel()
{
  dir = "MultiLayerModel";
  epsMethod       = multiLayerModel_EPS();
  isFinishMethod  = multiLayerModel_isFinish;
  needSizeMethod  = multiLayerModel_needSize;
  initModelMethod = multiLayerModel_init;  
  moveDirectoryMethod = multiLayerModel_moveDirectory;
}

static void morphoScaleModel()
{
  dir = "MorphoScaleModel";
  epsMethod = morphoScaleModel_EPS();
  isFinishMethod = morphoScaleModel_isFinish;
  moveDirectoryMethod = morphoScaleModel_moveDirectory;
  
  printf("not implemented morpho Model");
  exit(2);
}

static void concentricCircleModel()
{
  dir = "ConcentricCircleModel";
  makeDirectory(dir);
  moveDirectory(dir);
  epsMethod = concentricCircleModel_EPS();
  isFinishMethod      = concentricCircleModel_isFinish;
  moveDirectoryMethod = concentricCircleModel_moveDirectory;
}

bool models_isFinish()
{
  return (*isFinishMethod)();
}

//モデルを変更したときに,一度rootまで戻るので,再度ネストする用
void models_moveDirectory()
{
  moveDirectory(dir);
  (*moveDirectoryMethod)();
}

void models_setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_CYLINDER:
    circleModel();
    break;
  case LAYER:
    multiLayerModel();
    break;
  case MORPHO_SCALE:
    morphoScaleModel();
    break;
  case CONCENTRIC_CIRCLE:
    concentricCircleModel();
    break;
  }
  
  //ディレクトリの移動
//  models_moveDirectory();
}

double models_eps(double x, double y, enum MODE mode)
{
  double epsilon;
  switch(mode){
  case D_X :
    epsilon = (*epsMethod)(x, y, 1, 0);
  case D_Y :
    epsilon = (*epsMethod)(x, y, 0, 1);
  case D_XY :
    epsilon = (*epsMethod)(x, y, 1, 1);
  default:
    epsilon = (*epsMethod)(x, y, 1, 1);
  }
  
  return epsilon;
}

void models_needSize(int *x_nm,int *y_nm)
{
  (*needSizeMethod)(x_nm, y_nm);
}

void models_initModel()
{
  (*initModelMethod)();
}
