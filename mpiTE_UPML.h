#ifndef FDTD_TE_UPML_H
#define FDTD_TE_UPML_H
#include <complex.h>

extern void (* fdtdTE_upml_getUpdate(void))(void);
extern void (* fdtdTE_upml_getFinish(void))(void);
extern void (* fdtdTE_upml_getInit(void))(void);
extern double complex* fdtdTE_upml_getEx(void);
extern double complex* fdtdTE_upml_getEy(void);
extern double complex* fdtdTE_upml_getHz(void);

extern int fdtdTE_upml_getSubNx(void);
extern int fdtdTE_upml_getSubNy(void);
extern int fdtdTE_upml_getSubNpx(void);
extern int fdtdTE_upml_getSubNpy(void);
extern int fdtdTE_upml_getSubNcell(void);
extern void fdtdTE_upml_getSubFieldPositions(int *subNx,int *subNy,int *subNpx, int *subNpy);

extern double * fdtdTE_upml_getEps(void);
#endif
