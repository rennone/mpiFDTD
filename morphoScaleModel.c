#include "morphoScaleModel.h"
#include "field.h"
#include "function.h"
#include "parser.h"
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

double width_s[2];     //幅
double widthOnTopRate; //頂上に置ける幅の縮小率(基本0~1)
double thickness_s[2]; //厚さ
double n[2];
double ep_s[2];           //誘電率 = n*n*ep0_s
int layerNum;          //枚数
bool notsymmetry;      //左右比対称

static double eps(double x, double y, int col, int row)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double width = max(width_s[0], width_s[1]);
  double thick = thickness_s[0] + thickness_s[1];
  double height = thick*layerNum;  

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  int oy = fInfo_s.N_PY/2 - height/2;
  int ox = fInfo_s.N_PX/2;  
  double _x = x-ox;	//ox,oyを座標の原点に
  double _y = y-oy;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > (width/2+0.5) ||  _y < -0.5 || _y > height+0.5 )  
    return EPSILON_0_S;

  double s[2]={0,0}; //n1,n2それぞれの分割セルの数が入る

  for(double i=-16+0.5; i<16; i+=1){
    for(double j=-16+0.5; j<16; j+=1){
      double sx = _x + col*i/32.0; //細分化したセルの位置
      double sy = _y + row*j/32.0;

      //上下に飛び出ていないか確認
      if(sy < 0 || sy > height)
        continue;
      
      //thickで割ったあまり(double型なのでこんなやり方をしている)
      double modY = sy - floor(sy/thick)*thick;
      
      //境界上のときは両方の平均になる(普通は無い).
      if(modY == thickness_s[0]) {
        s[0] += 0.5*(fabs(sx) < width_s[0]/2);
        s[1] += 0.5*(fabs(sx) < width_s[1]/2);
        continue;
      }
      
      int k = (modY > thickness_s[0]); //どっちの屈折率にいるか調べる

      if (sx < 0 && notsymmetry)
        k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

      double p = 1-sy /height; //現在の高さを0~1で; 1=>頂上
      double _wid = width_s[k]*(p + (1-p)*widthOnTopRate);
      if(abs(sx) < _wid/2)
        s[k] +=1;
    }    
  }

  s[0] /= 32.0*32.0;
  s[1] /= 32.0*32.0;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep_s[0]*s[0] + ep_s[1]*s[1];
}

void readConfig()
{
  FILE *fp = NULL;
  if( !(fp = fopen("config.txt", "r")))
  {
    printf("cannot find config.txt of morphoScaleModel\n");
    exit(2);
  }

  int err;
  char buf[1024], tmp[1024];
  // 9文よみこみ
  for(int i=0; i<9; i++)
  {
    if( !parser_nextLine(fp, buf) )
    {
      printf("parse error, config.txt at MorphoScaleModel.c");
      exit(2);
    }
    
    if(i==6)
    {
      widthOnTopRate = strtod(buf, tmp);
    }
    else if(i == 7)
    {
      //最後の行はレイヤの数
      layerNum = atoi(buf);
    }
    else if(i==8)
    {
      notsymmetry = atoi(buf);
    }
    else
    {
      if( i%3 == 0)
        width_s[i/3] = field_toCellUnit(atoi(buf));
      else if( i%3 == 1)
        thickness_s[i/3] = field_toCellUnit(atoi(buf));
      else
      {
        n[i/3] = strtod(buf, tmp);
        ep_s[i/3] = n[i/3] * n[i/3] * EPSILON_0_S;
      }
    }
  }
  fclose(fp);
}

double ( *morphoScaleModel_EPS(void))(double, double, int, int)
{
  int rank, numProc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  if(rank == 0)
  {
    readConfig();
      
    for(int i=1; i<numProc; i++)
    {
      MPI_Send(&width_s, 2, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
      MPI_Send(&widthOnTopRate, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
      MPI_Send(&thickness_s, 2, MPI_DOUBLE, i, 2, MPI_COMM_WORLD);
      MPI_Send(&ep_s, 2, MPI_DOUBLE, i, 3, MPI_COMM_WORLD);
      MPI_Send(&layerNum, 1, MPI_INT, i, 4, MPI_COMM_WORLD);
      MPI_Send(&notsymmetry, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
    }
  } else {
    MPI_Status status;
    MPI_Recv(&width_s, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&widthOnTopRate, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&thickness_s, 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(&ep_s, 2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(&layerNum, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
    MPI_Recv(&notsymmetry, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
  }  
  
  char str[1024];
  sprintf(str, "w%d_%d_d%d_%d",
          (int)field_toPhisycalUnit(width_s[0]), (int)field_toPhisycalUnit(width_s[1]),
          (int)field_toPhisycalUnit(thickness_s[0]), (int)field_toPhisycalUnit(thickness_s[1])
    );  
  makeDirectory(str);
  moveDirectory(str);
  
  char str2[1024];
  sprintf(str2, "n%.2lf-%.2lf_m%d_r%.2lf", n[0],n[1],layerNum, widthOnTopRate);
  makeDirectory(str2);
  moveDirectory(str2);

  return eps;
}
