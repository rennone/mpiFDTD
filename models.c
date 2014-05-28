#include "field.h"
#include "models.h"
#include "noModel.h"
#include "circleModel.h"
#include "multiLayerModel.h"
#include "morphoScaleModel.h"
#include "bool.h"
//#include "shelf.h"
//#include "nonshelf.h"

static double (*epsMethod)(double, double, int, int);
static bool (*isFinishMethod)(void);

static void noModel(void)
{
  //no material
  moveDirectory("NoModel");
  epsMethod = noModel_EPS();
  isFinishMethod = noModel_isFinish();
}

static void circleModel(void)
{
  moveDirectory("MieCylinderModel");
  //cylinder material whitch radius = lambda, origin = center of field
  epsMethod = circleModel_EPS(N_PX*0.5, N_PY*0.5, field_getLambda());
  isFinishMethod = circleModel_isFinish();
}

static void multiLayerModel()
{
  moveDirectory("MultiLayerModel");
  epsMethod = multiLayerModel_EPS();
  isFinishMethod = multiLayerModel_isFinish();
}

static void morphoScaleModel()
{
  moveDirectory("MorphoScaleModel");
  epsMethod = morphoScaleModel_EPS();
  isFinishMethod = morphoScaleModel_isFinish();
}

bool models_isFinish()
{
  return (*isFinishMethod)();
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
  }
}

double models_eps(double x, double y, enum MODE mode){
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
