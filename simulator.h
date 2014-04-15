#ifndef _SIMULATOR_H
#define _SIMULATOR_H
#include<complex.h>
#include "bool.h"
#include "models.h"
#include "field.h"

enum SOLVER
{
  TM_UPML_2D,
  TE_UPML_2D,
};

extern void simulator_init(FieldInfo field_info, enum MODEL model, enum SOLVER solver);
//extern void simulator_init(int width, int height , double h_u, int pml, double lambda, int step, enum MODEL model, enum SOLVER solver);
extern void simulator_calc(void);
extern bool simulator_isFinish(void);
extern void simulator_finish(void);
extern double complex* simulator_getDrawingData();

extern void simulator_getSubFieldPositions(int *subNx,int *subNy,int *subNpx, int *subNpy);
extern double *simulator_getEps();

#endif
