#include "morphoScaleModel.h"
#include "field.h"
#include "function.h"
#include "parser.h"
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
//モルフォ蝶の鱗粉モデル.

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラの厚さ
#define ST_THICK_NM_0 90
#define EN_THICK_NM_0 150
#define DELTA_THICK_NM_0 30

//空気の部分の厚さ
#define ST_THICK_NM_1 90
#define EN_THICK_NM_1 150
#define DELTA_THICK_NM_1 30

//ラメラの枚数
#define LAYER_NUM 11

//互い違い
#define ASYMMETRY false

//中心に以下の幅で軸となる枝を入れる
#define ST_BRANCH_NM 0
#define EN_BRANCH_NM 0
#define DELTA_BRANCH_NM 10

//屈折率
#define N_0 1.56

//#define N_0 8.4179 //serikon

//先端における横幅の割合
#define ST_EDGE_RATE 1.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.1

//ラメラの先端を丸める曲率 (1で四角形のまま, 0.0で最もカーブする)
#define CURVE 1.0

//エッジの角度をランダムに傾ける
#define RANDOM_EDGE_ANGLE true
static double edge_randomeness[LAYER_NUM];

static int width_nm     = ST_WIDTH_NM;
static int thickness_nm[2] = {ST_THICK_NM_0, ST_THICK_NM_1};
static int layerNum = LAYER_NUM;     //枚数
static int branch_width_nm = ST_BRANCH_NM; //枝の幅

static double branch_width_s; //枝の幅
static double width_s;        //幅
static double thickness_s[2]; //厚さ

static double ep_s;        //誘電率 = n*n*ep0

static double edge_width_rate = ST_EDGE_RATE;

static double c0; //2次関数の比例定数

// 傾きがa, 切片がb1とb2 (b1<b2) の２つの平行な直線に, (x,y)が含まれるか判定
static bool in2Line(double a, double b1, double b2, double x, double y)
{
  double y1 = a*x + b1;
  double y2 = a*x + b2;

  return ( y1 < y && y < y2 );
}

/*
static double calc_width(double sx, double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
//  double new_wid = wid*(p + (1-p)*edge_width_rate);
  double new_wid = wid-branch_width_s*2 + (p + (1-p)*edge_width_rate)*branch_width_s;

  //ラメラの下を基準とした位置を求める
  double dh = k==0 ? modY : modY - thickness_s[0];
  double c  = k==0 ? c0 : c1;

//互い違いの場合はdhを再計算
  if(ASYMMETRY && sx < 0){
    dh = (k==1 ? modY : modY - thickness_s[1]);
  }

  //2次関数で横幅を計算
  return c*pow((dh-thickness_s[k]/2),2) + new_wid;
  }*/

static double eps(double x, double y, int col, int row)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  double width = width_s;
  double thick = thickness_s[0] + thickness_s[1];
  double height = thick*layerNum;  

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  int oy = fInfo_s.N_PY/2 - height/2;
  int ox = fInfo_s.N_PX/2;  
  double _x = x-ox;	//ox,oyを座標の原点に
  double _y = y-oy;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  double max_width = 0.5*(width*sqrt(2.0)/2.0 + thickness_s[0]); //45°の時に最も横幅が大きくなるのでそれを基準にする
    if( fabs(_x) > (max_width + 0.5) )
    return EPSILON_0_S;

  int num = -1;
  double a1, b1, b2;  
  for( int i=0; i<layerNum; i++){
    //エッジを回転させたときの, 二つの直線を求める.
    double edge_oy = thick*i + 0.5*thickness_s[0]; //エッジの中心位置
    double bm_b = edge_oy - 0.5*thickness_s[0]/cos(edge_randomeness[i]); //切片
    double tp_b = edge_oy + 0.5*thickness_s[0]/cos(edge_randomeness[i]);
    double a = _x < 0 ? -tan(edge_randomeness[i]) : tan(edge_randomeness[i]);

    // 1セルの正方形が重なるかを返す(4点のどれかが内側にあれば, 重なっていると判断. 2直線の幅は1セルよりも大きくするのでokなはず)
    if( in2Line(a, bm_b, tp_b, _x-0.5, _y-0.5) || in2Line(a, bm_b, tp_b, _x-0.5, _y+0.5) ||
        in2Line(a, bm_b, tp_b, _x+0.5, _y-0.5) || in2Line(a, bm_b, tp_b, _x+0.5, _y+0.5))
    {
      num = i;
      a1 = a;
      b1 = bm_b;
      b2 = tp_b;
      break;
    }
  }
  
  double s=0; //n1の分割セルの数が入る  
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/32.0; //細分化したセルの位置
      double sy = _y + row*j/32.0;

      double p = 1 - sy/height;
      //枝の部分にあるか
      if(abs(sx) < branch_width_s*(p + (1-p)*edge_width_rate)){
        s += 1;
        continue;
      }
      
      //どのエッジとも重ならなければ, 枝とだけ判定すればいい
      if(num <0)
        continue;

      //中心のxと符号が逆になる場合は, 傾きも逆になる.
      double a = sx*_x < 0 ? -a1 : a1;
      
      //TODO　エッジの幅を考慮してない
      if( in2Line(a, b1, b2, sx, sy) ){
        double sqrWidth = sx*sx*(1+a*a);
        double offset = (b2 - (sy-a*sx))*sin(edge_randomeness[num]);
        if( sqrWidth < pow(width/2 + offset, 2))
          s+=1;
      }
    }
  }
  s /= split*split;
  return EPSILON_0_S*(1-s) + ep_s*s;
}

double ( *morphoScaleModel_EPS(void))(double, double, int, int)
{
  return eps;
}

//正しいディレクトリまで移動.
void morphoScaleModel_moveDirectory()
{
  if(ASYMMETRY){
    makeDirectory("asymmetry");
    moveDirectory("asymmetry");
  } else {
    makeDirectory("symmetry");
    moveDirectory("symmetry");
  }  
  char buf[512];
  // make folder by index of reflaction 
  sprintf(buf,"n_%.2lf", N_0);
  makeDirectory(buf);
  moveDirectory(buf);

  sprintf(buf,"curve_%.2lf", CURVE);
  makeDirectory(buf);
  moveDirectory(buf);

  sprintf(buf, "thick%d_%d_layer%d_edge%.1lf_branch%d",
          thickness_nm[0], thickness_nm[1], layerNum, edge_width_rate, branch_width_nm);
  makeDirectory(buf);
  moveDirectory(buf);
}

void morphoScaleModel_init()
{
  width_s     = field_toCellUnit(width_nm);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[1]);
  ep_s = N_0*N_0*EPSILON_0_S;

  branch_width_s = field_toCellUnit(branch_width_nm);
  
  c0 = 4*width_s*(CURVE-1)/thickness_s[0]/thickness_s[0];

  for(int i=0; i<layerNum; i++)
    edge_randomeness[i] = 45.0*M_PI/180.0;
}

//構造を一つ進める
static bool nextStructure()
{
  thickness_nm[0] += DELTA_THICK_NM_0;
  thickness_nm[1] += DELTA_THICK_NM_1;

  if(thickness_nm[0] > EN_THICK_NM_0)
  {
    thickness_nm[0] = ST_THICK_NM_0;
    thickness_nm[1] = ST_THICK_NM_1;

    edge_width_rate += DELTA_EDGE_RATE;
    if(edge_width_rate > EN_EDGE_RATE)
    {
      edge_width_rate = ST_EDGE_RATE;

      branch_width_nm += DELTA_BRANCH_NM;
      if(branch_width_nm > EN_BRANCH_NM)
      	{
	  printf("there are no models which hasn't been simulated yet\n");     
	  return true;
	}
    }
  }
  return false;  
}

bool morphoScaleModel_isFinish(void)
{
  return nextStructure();
}

void morphoScaleModel_needSize(int *x, int *y)
{
  *x = width_nm + branch_width_nm;
  *y = (thickness_nm[0]+thickness_nm[1])*layerNum;
}

/*
static void readConfig()
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
    printf("%s",buf);
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
      asymmetry = atoi(buf);
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

  printf("========MorphoScaleModel=======\n");
  printf("width(%.1lf,%.1lf) thickness(%.1lf,%.1lf) n(%.3lf,%.3lf)\n",width_s[0], width_s[1], thickness_s[0], thickness_s[1], n[0], n[1]);
  printf("LayerNum=%d, rate=%.3lf\n",layerNum, widthOnTopRate);
  printf("==============================\n");
}
*/
/*
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
      MPI_Send(&asymmetry, 1, MPI_INT, i, 5, MPI_COMM_WORLD);
      MPI_Send(&n, 2, MPI_DOUBLE, i, 6, MPI_COMM_WORLD);
    }
  } else {
    MPI_Status status;
    MPI_Recv(&width_s, 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(&widthOnTopRate, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(&thickness_s, 2, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
    MPI_Recv(&ep_s, 2, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
    MPI_Recv(&layerNum, 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
    MPI_Recv(&asymmetry, 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
    MPI_Recv(&n, 2, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status);
  }  
  return eps;
}
*/
