#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
#include "concentricCircleModel.h"
#include "multiLayerModel.h"
#include "morphoScaleModel.h"
#include "bool.h"
//#include "shelf.h"
//#include "nonshelf.h"

static double (*epsMethod)(double, double, int, int);
static bool (*isFinishMethod)(void);
static void (*moveDirectoryMethod)(void);
static char *dir;
static void noModel(void)
{
  //no material
  moveDirectory("NoModel");
  epsMethod = noModel_EPS();
  isFinishMethod = noModel_isFinish;

}

static void circleModel(void)
{
  dir = "MieCylinderModel";
  moveDirectory(dir);
  //cylinder material whitch radius = lambda, origin = center of field
  epsMethod = circleModel_EPS(N_PX*0.5, N_PY*0.5, field_getLambda());  
  isFinishMethod = circleModel_isFinish;
  
  moveDirectoryMethod = circleModel_moveDirectory;
}

static void multiLayerModel()
{
  moveDirectory("MultiLayerModel");
  epsMethod = multiLayerModel_EPS();  
  isFinishMethod = multiLayerModel_isFinish;

  printf("not implemented MoveDirectoryMethod");
  exit(2);
}

static void morphoScaleModel()
{
  dir = "MorphoScaleModel";
  moveDirectory(dir);
  epsMethod = morphoScaleModel_EPS();
  isFinishMethod = morphoScaleModel_isFinish;
  moveDirectoryMethod = morphoScaleModel_moveDirectory;
}

static void concentricCircleModel()
{
  dir = "ConcentricCircleModel";
  makeDirectory(dir);
  moveDirectory(dir);
  epsMethod = concentricCircleModel_EPS();
//  epsMethod           = conc_eps;//concentricCircleModel_EPS();
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

void setModel(enum MODEL model)
{
  switch(model){
  case NO_MODEL:
    noModel();
    break;
  case MIE_CYLINDER:
    circleModel();
    break;
  case SHELF :
  case NONSHELF:
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
   (*moveDirectoryMethod)();
}

double models_eps(double x, double y, enum MODE mode)
{
  double epsilon = EPSILON_0_S;
  switch(mode){
  case D_X :
    epsilon = (*epsMethod)(x, y, 1, 0);
  case D_Y :
    epsilon = (*epsMethod)(x, y, 0, 1);
  case D_XY :
    epsilon = (*epsMethod)(x, y, 1, 1);
  }
  return epsilon;
}
