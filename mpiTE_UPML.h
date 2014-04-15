#ifndef FDTD_TE_UPML_H
#define FDTD_TE_UPML_H
#include <complex.h>

extern void (* mpi_fdtdTE_upml_getUpdate(void))(void);
extern void (* mpi_fdtdTE_upml_getFinish(void))(void);
extern void (* mpi_fdtdTE_upml_getInit(void))(void);
extern double complex* mpi_fdtdTE_upml_getEx(void);
extern double complex* mpi_fdtdTE_upml_getEy(void);
extern double complex* mpi_fdtdTE_upml_getHz(void);

extern int mpi_fdtdTE_upml_getSubNx(void);
extern int mpi_fdtdTE_upml_getSubNy(void);
extern int mpi_fdtdTE_upml_getSubNpx(void);
extern int mpi_fdtdTE_upml_getSubNpy(void);
extern int mpi_fdtdTE_upml_getSubNcell(void);
extern void mpi_fdtdTE_upml_getSubFieldPositions(int *subNx,int *subNy,int *subNpx, int *subNpy);

extern double * mpi_fdtdTE_upml_getEps(void);
#endif
