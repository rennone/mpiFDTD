#include "noModel.h"
#include "field.h"

//
static double eps(double x, double y, int col, int row)
{
  return EPSILON_0_S;
}

static void output(FILE *fp, double complex* data)
{
  printf("noModel doesn't output anythig data \n");
  fclose(fp);
}

double (*noModel_EPS(void))(double, double, int, int)
{
  return eps;
}

void (*noModel_output(void))(FILE *, double complex*)
{
  return output;
}

bool noModel_isFinish(void)
{
  return true;
}



