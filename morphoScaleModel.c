#include "morphoScaleModel.h"
#include "field.h"
#include "function.h"
#include "parser.h"
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
//モルフォ蝶の鱗粉モデル.

#define LEFT  false
#define RIGHT true

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
#define ST_LAYER_NUM 11
#define EN_LAYER_NUM 11
#define DELTA_LAYER_NUM 1

//互い違い
#define ASYMMETRY true

//中心に以下の幅で軸となる枝を入れる
#define ST_BRANCH_NM 50
#define EN_BRANCH_NM 50
#define DELTA_BRANCH_NM 10

//屈折率
#define N_0 1.56

//#define N_0 8.4179 //serikon

//先端における横幅の割合
#define ST_EDGE_RATE 0.5
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.1

//ラメラの先端を丸める曲率 (1で四角形のまま, 0.0で最もカーブする)
#define CURVE 1.0

//エッジの角度をランダムに傾ける
#define RANDOM_EDGE_ANGLE true
static double edge_randomeness[2][EN_LAYER_NUM];

static int width_nm     = ST_WIDTH_NM;
static int thickness_nm[2] = {ST_THICK_NM_0, ST_THICK_NM_1};
static int layerNum = ST_LAYER_NUM;     //枚数
static int branch_width_nm = ST_BRANCH_NM; //枝の幅

static double branch_width_s; //枝の幅
static double width_s;        //幅
static double thickness_s[2]; //厚さ
static double height_s; //全体の高さ.
static double ox_s, oy_s; //原点

static double ep_s;        //誘電率 = n*n*ep0

static double edge_width_rate = ST_EDGE_RATE;

static double c0; //2次関数の比例定数

//1ラメラの式
typedef struct Lamela
{
  double a, b1, b2; //傾き, 下の切片, 上の切片
  int id; //ラメラのid
} Lamela;

static bool inLamela( Lamela *lamera, double x, double y)
{
  return ( lamera->a*x + lamera->b1 <= y && y <= lamera->a*x + lamera->b2 );
}

/*
static double calc_width(double sx, double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
  double new_wid = wid-branch_width_s*2 + (p + (1-p)*edge_width_rate)*branch_width_s;

  //ラメラの下を基準とした位置を求める
  double dh = k==0 ? modY : modY - thickness_s[0];

//互い違いの場合はdhを再計算
  if(ASYMMETRY && sx < 0){
    dh = (k==1 ? modY : modY - thickness_s[1]);
  }

  //2次関数で横幅を計算
  return c*pow((dh-thickness_s[k]/2),2) + new_wid;
  }*/

//
static double calcLamela(double lft, double rht, double btm, double top, Lamela *lam)
{
  lam->id = -1;
  bool signX = lft<0 ? LEFT : RIGHT;
  bool reverse = (signX==LEFT) && ASYMMETRY ? true : false;
  double thick = thickness_s[0]+thickness_s[1];
  double half_thick = reverse ? 0.5*thickness_s[1] : 0.5*thickness_s[0];
  double offset_y = reverse ? thickness_s[0] + 0.5*thickness_s[1] : 0.5*thickness_s[0];
  for(int i=0; i<layerNum; i++)
  {
    double oy = thick*i + offset_y ; //エッジの中心位置;

    lam->a = tan(edge_randomeness[signX][i]); //傾き
    lam->b1 = oy - half_thick/cos(edge_randomeness[signX][i]); //切片
    lam->b2 = oy + half_thick/cos(edge_randomeness[signX][i]);

    // 1セルの正方形が重なるかを返す(4点のどれかが内側にあれば, 重なっていると判断. 2直線の幅は1セルよりも大きくするのでokなはず)
    if( inLamela(lam, lft, btm) || inLamela(lam, rht, btm) ||
        inLamela(lam, lft, top) || inLamela(lam, rht, top))
    {
      lam->id = i;
      break;
    }
  }
}

static double eps(double x, double y, int col, int row)
{
  double thick = thickness_s[0] + thickness_s[1];
  double _x = x-ox_s;	//ox,oyを座標の原点に
  double _y = y-oy_s;

  //ラメラの中心で回転するので,本来の横幅よりも長くなる可能性がある.
  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > 0.5*sqrt(width_s*width_s + thickness_s[0]*thickness_s[0])+0.5 )
    return EPSILON_0_S;

  Lamela lams[2];

  if( fabs(_x) > 0)
  {
    calcLamela(_x-0.5, _x+0.5, _y-0.5, _y+0.5, _x<0 ? &lams[LEFT] : &lams[RIGHT]);
  } else
  {
    calcLamela(_x-0.5,      0, _y-0.5, _y+0.5, &lams[LEFT]);
    calcLamela( 0    , _x+0.5, _y-0.5, _y+0.5, &lams[RIGHT]);
  }

  double s=0; //n1の分割セルの数が入る  
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/split; //細分化したセルの位置
      double sy = _y + row*j/split;
      double p = 1 - sy/height_s;
      double width_p = (p + (1-p)*edge_width_rate);
      //枝の部分にあるか
      if(fabs(sx) < branch_width_s*width_p && sy>=0 && sy <= height_s){
        s += 1;
        continue;
      }
      
      bool signX = sx < 0 ? LEFT : RIGHT;
      Lamela *lam = &lams[signX];
      
      //どのエッジとも重ならなければ, 枝とだけ判定すればいい
      if( lam->id < 0)
        continue;

      //ラメラのエッジの変化は, 切片に写像して考える.      
      p = 1 - (sy-lam->a*sx)/height_s;
      width_p = (p + (1-p)*edge_width_rate);
      
      if( inLamela(lam, sx, sy) )
      {
        double sqrWidth = pow(sx,2)*(1+pow(lam->a,2) );
        double rad = edge_randomeness[signX][lam->id];        
        double offset = (signX == LEFT ? -1 : 1) *( 0.5*(lam->b2+lam->b1) - (sy-lam->a*sx) ) * sin( rad );
        if( sqrWidth < pow( width_p*(width_s/2 + offset), 2))
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
  branch_width_s = field_toCellUnit(branch_width_nm);
  
  height_s = (thickness_s[0] + thickness_s[1])*layerNum;  
  ep_s = N_0*N_0*EPSILON_0_S;
  
  c0 = 4*width_s*(CURVE-1)/thickness_s[0]/thickness_s[0];

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  oy_s = fInfo_s.N_PY/2 - height_s/2;
  ox_s = fInfo_s.N_PX/2;

  double TO_RAD = M_PI/180.0;
  for(int i=0; i<layerNum; i++)
  {
    edge_randomeness[RIGHT][i] = (rand()%20-10)*TO_RAD;
    edge_randomeness[LEFT][i]  = (rand()%20-10)*TO_RAD;
  }
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
