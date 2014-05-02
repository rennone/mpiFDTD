#include "function.h"

double dbilinear(double *p, double x, double y, int width, int height)
{
  int i = floor(x);
  int j = floor(y);
  double dx = x - i;
  double dy = y - j;
  int index = i*height + j;
  return p[index]*(1.0-dx)*(1.0-dy)
       + p[index+height]*dx*(1.0-dy)
       + p[index+1]*(1.0-dx)*dy
       + p[index+height+1]*dx*dy;
}

FILE* openFile(const char* file_name)
{
  FILE *fp;
  if( (fp=fopen(file_name, "w") ) == NULL )
  {
    printf("cannot open file %s \n", file_name);
    exit(2);
  }
  return fp;
}
