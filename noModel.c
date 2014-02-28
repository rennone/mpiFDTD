#include "noModel.h"
#include "field.h"

static double eps(double x, double y, int col, int row)
{
  return EPSILON_0_S;
}

double (*noModel_EPS(void))(double, double, int, int)
{
  return eps;
}

static void output(FILE *fp, double complex* data)
{
  printf("noModel doesn't output anythig data \n");
  fclose(fp);
}


void (*noModel_output(void))(FILE *, double complex*)
{
  return output;
}



