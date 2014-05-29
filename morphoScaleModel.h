#ifndef MORPHO_SCALE_H
#define MORPHO_SCALE_H

#include "bool.h"
extern void morphoScaleModel_init();
extern double ( *morphoScaleModel_EPS(void))(double, double, int, int);
extern bool morphoScaleModel_isFinish(void);

extern void morphoScaleModel_moveDirectory(void);
#endif
