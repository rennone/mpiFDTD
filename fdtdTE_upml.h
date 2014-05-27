#ifndef FDTD_TE_UPML_H
#define FDTD_TE_UPML_H
#include "myComplex.h"

extern void (* fdtdTE_upml_getUpdate(void))(void);
extern void (* fdtdTE_upml_getFinish(void))(void);
extern void (* fdtdTE_upml_getInit(void))(void);
extern void (* fdtdTE_upml_getReset(void))(void);
extern double complex* fdtdTE_upml_getEx(void);
extern double complex* fdtdTE_upml_getEy(void);
extern double complex* fdtdTE_upml_getHz(void);


extern double* fdtdTE_upml_getEps();
#endif
