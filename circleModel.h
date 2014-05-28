#ifndef _CIRCLE_H
#define _CIRCLE_H

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include "bool.h"
extern double (*circleModel_EPS(double x, double y, double radius))(double, double, int,int);
extern void (*circleModel_output(void))(FILE *, double complex*);

extern bool circleModel_isFinish(void);
#endif
