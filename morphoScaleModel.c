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
#define ST_THICK_NM_0 80
#define EN_THICK_NM_0 170
#define DELTA_THICK_NM_0 30

//空気の部分の厚さ
#define ST_THICK_NM_1 80
#define EN_THICK_NM_1 170
#define DELTA_THICK_NM_1 30

//ラメラの枚数
#define ST_LAYER_NUM 4
#define EN_LAYER_NUM 12
#define DELTA_LAYER_NUM 2

//互い違い
#define ASYMMETRY true

//中心に以下の幅で軸となる枝を入れる
#define ST_BRANCH_NM 0
#define EN_BRANCH_NM 50
#define DELTA_BRANCH_NM 10

//屈折率
#define N_0 1.56

//#define N_0 8.4179 //serikon

//先端における横幅の割合
#define ST_EDGE_RATE 0.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.5

//ラメラの先端を丸める曲率 (0で四角形のまま, 1.0で最もカーブする)
#define CURVE 0.2

//エッジの角度をランダムに傾ける
#define RANDOMNESS 20
#define RANDOM_SEED 0 //各プロセスで同じ角度になるようにseedを固定する
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

static double calcWidth(double sx, double sy, Lamela *lam)
{
  //x座標により, width_rateが変わるのを防ぐため.
  //ラメラのエッジの変化は, (sx,sy)を切片に写像して考える.
  bool signX = sx < 0 ? LEFT : RIGHT;
  double rad = edge_randomeness[signX][lam->id];        
  double offset = (signX == LEFT ? -1 : 1) *( 0.5*(lam->b2+lam->b1) - (sy-lam->a*sx) ) * sin( rad );

  double p = 1 - (sy-lam->a*sx)/height_s;
  double width_p = (p + (1-p)*edge_width_rate);
  double b = (sy-sx*lam->a);
  double x  =  (b-lam->b1)/(lam->b2-lam->b1) - 0.5;

  double half_width = width_s/2.0;
  // 両端(x=+-0.5)の時に, width_s/2.0*CURVE 長さが減るように係数を設定
  double c0 = 4.0*half_width*CURVE;
  
  //丸まったところ以外だけ長さが短くなるようにする.
  // p*(丸みをのぞいた横幅) + 丸みまった部分 + 傾きによるオフセット
  return width_p*( half_width*(1-CURVE) ) + c0*(0.25 - x*x) + offset;
}

static void calcLamela(double lft, double rht, double btm, double top, Lamela *lam)
{
  lam->id = -1;
  bool signX = lft<0 ? LEFT : RIGHT;
  bool reverse = (signX==LEFT) && ASYMMETRY ? true : false;
  double thick = thickness_s[0]+thickness_s[1];
  double half_thick = 0.5*thickness_s[0];//reverse ? 0.5*thickness_s[1] : 0.5*thickness_s[0];
  double offset_y = reverse ? thickness_s[1] + half_thick : half_thick;
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
  double _x = x-ox_s;	//ox,oyを座標の原点に
  double _y = y-oy_s;

  //ラメラの中心で回転するので,本来の横幅よりも長くなる可能性がある.
  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > 0.5*sqrt( pow(width_s,2)+pow(thickness_s[0],2) ) + 0.5 )
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
      
      //枝の部分にあるか判定
      if(sy>=0 && sy <= height_s){
        //ブランチはラメラとは無関係なのでpはそのまま計算.
        double p = 1 - sy/height_s;
        double width_p = (p + (1-p)*edge_width_rate);
        if(fabs(sx) < branch_width_s/2.0*width_p){
          s += 1;
          continue;
        }
      }
      
      bool signX = sx < 0 ? LEFT : RIGHT;
      Lamela *lam = &lams[signX];
      
      //どのエッジとも重ならなければ, 枝とだけ判定すればいい
      if( lam->id < 0)
        continue;
      
      if( inLamela(lam, sx, sy) )
      {
        double sqrWidth = pow(sx,2)*(1+pow(lam->a,2) );
        double newWidth = calcWidth(sx, sy, lam);
        if( sqrWidth < pow( newWidth, 2))
          s+=1;
      }
    }
  }
  s /= split*split;
  return EPSILON_0_S*(1-s) + ep_s*s;
}

//空気の層がエッジよりも厚みが少ないとき, 2つのラメラが重なる可能性があるのでチェック
static bool cross2Lamela(double wid, double thick1, double thick2, double rad)
{
  double length = sin(rad/2.0)*wid + (cos(rad/2.0)-1.0)*thick1;
  //上のラメラが-RANDOMNESS/2 , 下のラメラが+RANDOMNESS/2の時に
  //飛び出た分が空気の層の厚みを超えなければok
  return length > thick2;
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

  sprintf(buf, "randome%d_%d", RANDOMNESS, RANDOM_SEED );
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

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  oy_s = fInfo_s.N_PY/2 - height_s/2;
  ox_s = fInfo_s.N_PX/2;

  srand(RANDOM_SEED);
  double TO_RAD = M_PI/180.0;
  for(int i=0; i<layerNum; i++)
  {
    edge_randomeness[RIGHT][i] = (rand()%(RANDOMNESS+1)-RANDOMNESS/2)*TO_RAD;
    edge_randomeness[LEFT][i]  = (rand()%(RANDOMNESS+1)-RANDOMNESS/2)*TO_RAD;
  }
}

double ( *morphoScaleModel_EPS(void))(double, double, int, int)
{
  
  /*
  if(cross2Lamela())
  {
    printf("skip this case because perhaps two lamela will be cross \n");
    bool res = nextStructure();
    while ( !res && cross2Lamela() )
    {
      printf("skip this case because perhaps two lamela will be cross \n");
      res = nextStructure();
    }
    if(res)
    {
      printf("can not try this simulation because all parameter has possible that two lamela will be cross \n");
      exit(0);
    }
    }*/
  if(cross2Lamela(EN_WIDTH_NM, ST_THICK_NM_0, ST_THICK_NM_1, RANDOMNESS*M_PI/180.0))
  {
    printf("can not try this simulation because there are possible that two lamela will be cross \n");
    exit(0);
  }
  return eps;
}

bool morphoScaleModel_isFinish(void)
{
  bool res = nextStructure();
  /*
  while ( !res && cross2Lamela() )
  {
    printf("skip this case because perhaps two lamela will be cross \n");
    res = nextStructure();
    }*/
  return res;
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
