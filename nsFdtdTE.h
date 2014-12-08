#ifndef _NS_FDTD_TE_H
#define _NS_FDTD_TE_H
#include "myComplex.h"
#include "solver.h"

extern void(*nsFdtdTE_getUpdate(void))(void);
extern void(*nsFdtdTE_getFinish(void))(void);
extern void(*nsFdtdTE_getReset(void))(void);
extern void(* nsFdtdTE_getInit(void))(void);

extern dcomplex* nsFdtdTE_getHzx(void);
extern dcomplex* nsFdtdTE_getHzy(void);
extern dcomplex* nsFdtdTE_getHz(void);
extern dcomplex* nsFdtdTE_getEy(void);
extern dcomplex* nsFdtdTE_getEx(void);

extern double* nsFdtdTE_getEps();

extern Solver* nsFdtdTE_getSolver(void);
#endif
