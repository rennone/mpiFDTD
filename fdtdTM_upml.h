#ifndef FDTD_TM_UPML_H
#define FDTD_TM_UPML_H
#include "myComplex.h"

extern void (* fdtdTM_upml_getUpdate(void))(void);
extern void (* fdtdTM_upml_getFinish(void))(void);
extern void (* fdtdTM_upml_getInit(void))(void);
extern dcomplex* fdtdTM_upml_getHx(void);
extern dcomplex* fdtdTM_upml_getHy(void);
extern dcomplex* fdtdTM_upml_getEz(void);

extern double* fdtdTM_upml_getEps();
#endif
