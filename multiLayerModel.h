#ifndef MULTI_LAYER_MODEL_H
#define MULTI_LAYER_MODEL_H

#include "bool.h"
#include <stdio.h>
#include <stdlib.h>
extern double ( *multiLayerModel_EPS(void))(double, double, int, int);
extern bool multiLayerModel_isFinish(void);
//extern void (*multiLayerModel_output(void))(FILE *, double complex*);

#endif
