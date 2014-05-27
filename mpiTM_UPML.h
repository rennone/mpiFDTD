#ifndef MPI_FDTD_TM_UPML_H
#define MPI_FDTD_TM_UPML_H
#include <complex.h>

extern void (* mpi_fdtdTM_upml_getUpdate(void))(void);
extern void (* mpi_fdtdTM_upml_getFinish(void))(void);
extern void (* mpi_fdtdTM_upml_getReset(void))(void);
extern void (* mpi_fdtdTM_upml_getInit(void))(void);
extern double complex* mpi_fdtdTM_upml_getHx(void);
extern double complex* mpi_fdtdTM_upml_getHy(void);
extern double complex* mpi_fdtdTM_upml_getEz(void);
extern double * mpi_fdtdTM_upml_getEps(void);

extern int mpi_fdtdTM_upml_getSubNx(void);
extern int mpi_fdtdTM_upml_getSubNy(void);
extern int mpi_fdtdTM_upml_getSubNpx(void);
extern int mpi_fdtdTM_upml_getSubNpy(void);
extern int mpi_fdtdTM_upml_getSubNcell(void);

#endif
