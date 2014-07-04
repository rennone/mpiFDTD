#include <math.h>
#include "circleModel.h"
#include "field.h"
#include "function.h"

static double radius;
static double epsilon;
static double posx;
static double posy;

static double eps(double, double, int, int);

double (*circleModel_EPS(double x, double y, double r))(double, double, int , int)
{
  radius = field_toCellUnit(500);
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  posx   = fInfo_s.N_PX/2;
  posy   = fInfo_s.N_PY/2;
  
  epsilon = 1.6*1.6*EPSILON_0_S;
  
  return eps;
}

//col : D_Xモード, row : D_Yモード
static double eps(double x, double y, int col, int row)
{
  if(x < N_PML || y < N_PML || x > N_X+N_PML || y > N_Y + N_PML)
    return EPSILON_0_S;

  double dx = x-posx;
  double dy = y-posy;
  //2乗距離
  double len = dx*dx+dy*dy;

  //中心との距離がr+1セル以上なら,そのセルは完全に媒質の外 
  if(len >= (radius+1)*(radius+1))
    return EPSILON_0_S;

  //中心との距離がr-1セル以下なら,そのセルは完全に媒質の外 
  if(len <= (radius-1)*(radius-1))
    return epsilon;

  //さらに32*32分割し媒質内と媒質外の数を求めepsilonを決定する
  double split = 32;
  double half_split = split/2;
  double sum=0;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      if(pow(dx+col*i/split, 2.0) + pow(dy+row*j/split, 2.0) <= radius*radius)
	sum+=1;
    }
  }
  
  sum /= split*split;
  return epsilon*sum + EPSILON_0_S*(1-sum);
}


static void output(FILE *fp, double complex* data)
{
  //半径の1.2倍の位置のデータを保存
  double observation = 1.2*radius;
  int ang;
  for(ang=0; ang <360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double x = observation*cos(rad)+N_PX/2.0;
    double y = observation*sin(rad)+N_PY/2.0;
    double norm = cnorm(cbilinear(data,x,y,N_PX,N_PY));
    fprintf(fp, "%d %lf \n", 180-ang, norm);   
  }  
  fclose(fp);  
  printf("output end \n");
}

void (*circleModel_output(void))(FILE *, double complex*)
{
  return output;
}

bool circleModel_isFinish()
{
  radius += field_toCellUnit(200);
  return radius > field_toCellUnit(2000);
}

void circleModel_moveDirectory()
{
  
}
