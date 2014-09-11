#ifndef _FDTD_TM_H
#define _FDTD_TM_H
#include "myComplex.h"
#include "solver.h"

extern Solver* fdtdTM_getSolver(void);

extern void(*fdtdTM_getUpdate(void))(void);
extern void(*fdtdTM_getFinish(void))(void);
extern void(*fdtdTM_getReset(void))(void);
extern void(* fdtdTM_getInit(void))(void);

extern dcomplex* fdtdTM_getEzx(void);
extern dcomplex* fdtdTM_getEzy(void);
extern dcomplex* fdtdTM_getEz(void);
extern dcomplex* fdtdTM_getHy(void);
extern dcomplex* fdtdTM_getHx(void);

extern double* fdtdTM_getEps();
#endif
