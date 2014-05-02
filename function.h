#ifndef _FUNCTION_H
#define _FUNCTION_H
#include "myComplex.h"
#include <stdio.h>

#define max(a,b) (a > b ? a : b)
#define min(a,b) (a < b ? a : b)

extern double dbilinear(double *p, double x, double y, int width, int height);

extern FILE* openFile(const char* file_name);

// for(i=1..N_PX-1)
//   for(j=1..N_PY-1) と同じ
#define FAST_FOR_FOR(k, fInfo_s) \
  for(int k=fInfo_s.N_PY+1, last = fInfo_s.N_CELL-fInfo_s.N_PY; k<last; k+=2) \
    for(int endRow = k+fInfo_s.N_PY-2; k<endRow; k++)


#endif
