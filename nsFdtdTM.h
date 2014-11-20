#ifndef _NS_FDTD_TM_H
#define _NS_FDTD_TM_H
#include "myComplex.h"
#include "solver.h"

extern void(*nsFdtdTM_getUpdate(void))(void);
extern void(*nsFdtdTM_getFinish(void))(void);
extern void(*nsFdtdTM_getReset(void))(void);
extern void(* nsFdtdTM_getInit(void))(void);

extern dcomplex* nsFdtdTM_getEzx(void);
extern dcomplex* nsFdtdTM_getEzy(void);
extern dcomplex* nsFdtdTM_getEz(void);
extern dcomplex* nsFdtdTM_getHy(void);
extern dcomplex* nsFdtdTM_getHx(void);

extern double* nsFdtdTM_getEps();

extern Solver* nsFdtdTM_getSolver(void);
#endif
