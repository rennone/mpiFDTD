#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

//屈折率
#define N_0 1.0
#define N_1 1.56
//serikon
//#define N_1 8.4179

//異方性を入れるかのフラグ
#define UNIAXIAL false
#define N_0_X 1.0
#define N_1_X 1.1

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラ1の厚さ
#define ST_THICK_NM_1 70
#define EN_THICK_NM_1 130
#define DELTA_THICK_NM_1 30

//ラメラ0の厚さ
#define ST_THICK_NM_0 (ST_THICK_NM_1 + 10)
#define EN_THICK_NM_0 (ST_THICK_NM_1 + 50)
#define DELTA_THICK_NM_0 10

//ラメラの枚数
#define ST_LAYER_NUM 4
#define EN_LAYER_NUM 4
#define DELTA_LAYER_NUM 2

//互い違い => 左右で n_0 n_1を入れ替え
#define ASYMMETRY false

//USE_GAP flag 
//左右でずらす => DELTA_LEFT_GAP_Y ~ thickness_nm まで変化する.
#define USE_GAP false
#define DELTA_LEFT_GAP_Y 10

//中心に以下の幅で軸となる枝を入れる => 軸の屈折率はN_1になる
#define ST_BRANCH_NM 0
#define EN_BRANCH_NM 0
#define DELTA_BRANCH_NM 50

//先端における横幅の割合
#define ST_EDGE_RATE 0.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 1.0

//ラメラの先端を丸める曲率 (0で四角形のまま, 1.0で最もカーブする)
#define CURVE 0.0

static int width_nm[2]     = {ST_WIDTH_NM, ST_WIDTH_NM};
static double width_s[2];     //幅

static int thickness_nm[2] = {ST_THICK_NM_0, ST_THICK_NM_1};
static double thickness_s[2]; //厚さ

static int branch_width_nm = ST_BRANCH_NM; //枝の幅
static double branch_width_s; //枝の幅

static int layerNum = ST_LAYER_NUM;     //枚数

//N_0 N_1から計算
static double ep_s[2];        //誘電率 = n*n*ep0
static double ep_x_s[2];      //異方性用のx方向の誘電率

static double edge_width_rate = ST_EDGE_RATE;

//DELTA_LEFT_GAP_Y
#if USE_GAP
static int left_gap_y_nm = DELTA_LEFT_GAP_Y;
#else
static int left_gap_y_nm = 0;
#endif

static double left_gap_y_s;

//CURVE から計算
static double c0, c1; //2次関数の比例定数

static double calc_width(double sx, double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
  double new_wid = (wid+branch_width_s)*(p + (1-p)*edge_width_rate);

  //ラメラの下を基準とした位置を求める
  double dh = k==0 ? modY : modY - thickness_s[0];
  double c  = k==0 ? c0 : c1;

//互い違いの場合はdhを再計算
  if(ASYMMETRY && sx < 0){
    dh = (k==1 ? modY : modY - thickness_s[1]);
  }

  //2次関数で横幅を計算
  return c*pow((dh-thickness_s[k]/2),2) + new_wid;
}

//col : D_Xモード row : D_Yモード
//x,yを中心に, 計算領域のセルと同じ大きさの領域を調べる
static double eps(double x, double y, int col, int row)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();

  double width = max(width_s[0], width_s[1]);
  double thick = thickness_s[0] + thickness_s[1];
  double height = thick*layerNum + left_gap_y_s; //オフセットがあれば行う

  //領域の中心から, 下にheight/2ずれた位置がレイヤの下部
  int oy = fInfo_s.N_PY/2 - height/2;
  int ox = fInfo_s.N_PX/2;
  double _x = x-ox;	//ox,oyを座標の原点に
  double _y = y-oy;

  //上下左右に飛び出ていないか確認(細分化したセルがあるため, 0.5の余白をとっている)
  if( fabs(_x) > (width/2+0.5) ||  _y < -0.5 || _y > height+0.5 )  
    return EPSILON_0_S;

  height = thick*layerNum; //元にもどす.
  double s[2]={0,0}; //n1,n2それぞれの分割セルの数が入る
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/split; //細分化したセルの位置
      double sy = _y + row*j/split;

      //左側はleft_gapだけずらす.
      if(sx < 0 && USE_GAP)
        sy -= left_gap_y_s;
      
      //上下に飛び出ていないか確認
      if(sy < 0 || sy > height)
        continue;

      double p = 1 - sy/height;
      //枝の部分
      if(fabs(sx) < branch_width_s*(p + (1-p)*edge_width_rate))
      {
        s[1] += 1;
        continue;
      }
      
      //thickで割ったあまり(double型なのでこんなやり方をしている)
      double modY = sy - floor(sy/thick)*thick;

      //境界上のときは両方の平均になる(普通は無い).
      if( modY == thickness_s[0]) {
        s[0] += 0.5*(fabs(sx) < width_s[0]/2);
        s[1] += 0.5*(fabs(sx) < width_s[1]/2);
        continue;
      }

      //どっちの屈折率にいるか調べる
      int k;
      if (sx < 0 && ASYMMETRY) {
        k = (modY < thickness_s[1]); //互い違いかつ左側は, 1が下にある
      } else {
        k = (modY > thickness_s[0]); //それ以外は0が下にある
      }
      
      double wid = calc_width(sx, sy, width_s[k], height, modY, k);
      
      if(fabs(sx) < wid/2)
        s[k] +=1;
    }    
  }

  s[0] /= split*split;
  s[1] /= split*split;
  if(UNIAXIAL && col == 0 && row == 1)
  {
    return EPSILON_0_S*(1-s[0]-s[1]) + ep_x_s[0]*s[0] + ep_x_s[1]*s[1];
  }
  else
  {
    return EPSILON_0_S*(1-s[0]-s[1]) + ep_s[0]*s[0] + ep_s[1]*s[1];
  }
}

double ( *multiLayerModel_EPS(void))(double, double, int, int)
{
  return eps;
}

/*
//構造を一つ進める
static bool nextStructure()
{
  left_gap_y_nm += DELTA_LEFT_GAP_Y;
  if( left_gap_y_nm >= (thickness_nm[0]+thickness_nm[1]) || !USE_GAP)
  {
    left_gap_y_nm = USE_GAP ? DELTA_LEFT_GAP_Y : 0;
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
	  branch_width_nm = ST_BRANCH_NM;
	  layerNum += DELTA_LAYER_NUM;
	  if( layerNum > EN_LAYER_NUM)
          {
            printf("there are no models which hasn't been simulated yet\n");
            return true;
          }
	}
      }
    }
  }
  return false;  
}
*/

//構造を一つ進める
static bool nextStructure()
{
  left_gap_y_nm += DELTA_LEFT_GAP_Y;
  if( left_gap_y_nm >= (thickness_nm[0]+thickness_nm[1]) || !USE_GAP){
  left_gap_y_nm = USE_GAP ? DELTA_LEFT_GAP_Y : 0;
    thickness_nm[0] += DELTA_THICK_NM_0;

    if(thickness_nm[0] > EN_THICK_NM_0){
      thickness_nm[1] += DELTA_THICK_NM_1;
      thickness_nm[0] = ST_THICK_NM_0;
    
      if(thickness_nm[1] > EN_THICK_NM_1){ 
	thickness_nm[1] = ST_THICK_NM_1;
	edge_width_rate += DELTA_EDGE_RATE;
      
	if(edge_width_rate > EN_EDGE_RATE) {
	  edge_width_rate = ST_EDGE_RATE;
	  branch_width_nm += DELTA_BRANCH_NM;
	
	  if(branch_width_nm > EN_BRANCH_NM) {
	    branch_width_nm = ST_BRANCH_NM;
	    layerNum += DELTA_LAYER_NUM;
	  
	    if( layerNum > EN_LAYER_NUM) {
	      printf("there are no models which hasn't been simulated yet\n");
	      return true;
	    }
	  }
	}
      }
    }
  }
  return false;  
}


bool multiLayerModel_isFinish(void)
{
  return nextStructure();
}

void multiLayerModel_needSize(int *x_nm, int *y_nm)
{
  (*x_nm) = max( width_nm[0], width_nm[1]) + branch_width_nm;

  //最後の項はgapの分(これは固定にしないと,gapによりフィールドの領域が変わるので図が変に見える).
  (*y_nm) = (thickness_nm[0]+thickness_nm[1])*layerNum + (USE_GAP ? thickness_nm[0]+thickness_nm[0] : 0) ;
}

void multiLayerModel_moveDirectory()
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

  if(UNIAXIAL){
    sprintf(buf,"uniaxial_n0y%.2lf_n0x%.2lf_n1x%.2lf_n1x%.2lf", N_0, N_0_X, N_1, N_1_X);
    makeDirectory(buf);
    moveDirectory(buf);
  }
  else
  {
    sprintf(buf,"n_%.2lf_%.2lf", N_0, N_1);
    makeDirectory(buf);
    moveDirectory(buf);
  }
  
  sprintf(buf,"curve_%.2lf", CURVE);
  makeAndMoveDirectory(buf);

  sprintf(buf, "width%d_%d", width_nm[0], width_nm[0]);
  makeAndMoveDirectory(buf);

  sprintf(buf, "thick%d_%d", thickness_nm[0], thickness_nm[1]);
  makeAndMoveDirectory(buf);
  
  sprintf(buf, "gap%d", left_gap_y_nm);
  makeAndMoveDirectory(buf);

  sprintf(buf, "layer%d", layerNum);
  makeAndMoveDirectory(buf);

  sprintf(buf, "edge%.1lf", edge_width_rate);
  makeAndMoveDirectory(buf);

  sprintf(buf, "branch%d", branch_width_nm);
  makeAndMoveDirectory(buf);
}

void multiLayerModel_init()
{
  width_s[0]     = field_toCellUnit(width_nm[0]);
  width_s[1]     = field_toCellUnit(width_nm[1]);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[1]);
  left_gap_y_s   = field_toCellUnit(left_gap_y_nm);
  
  ep_s[0] = N_0*N_0*EPSILON_0_S;
  ep_s[1] = N_1*N_1*EPSILON_0_S;
  ep_x_s[0] = N_0_X*N_0_X*EPSILON_0_S;
  ep_x_s[1] = N_1_X*N_1_X*EPSILON_0_S;
  
  branch_width_s = field_toCellUnit(branch_width_nm);
  
  c0 = -4*width_s[0]*CURVE/thickness_s[0]/thickness_s[0];
  c1 = -4*width_s[1]*CURVE/thickness_s[1]/thickness_s[1];
}
