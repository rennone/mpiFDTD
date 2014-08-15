#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

/*
ASYMMETRYがtrueの場合, ラメラ1,2が同じ幅じゃないと, 奇麗にに互い違いにならない(左右で幅が異なってしまうから).
 */
//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラの厚さ
#define ST_THICK_NM_0 60
#define EN_THICK_NM_0 160
#define DELTA_THICK_NM_0 10

#define ST_THICK_NM_1 60
#define EN_THICK_NM_1 160
#define DELTA_THICK_NM_1 10

//ラメラの枚数
#define ST_LAYER_NUM 4
#define EN_LAYER_NUM 8
#define DELTA_LAYER_NUM 1
//#define LAYER_NUM 4

//互い違い
#define ASYMMETRY true

//中心に以下の幅で軸となる枝を入れる => 軸の屈折率はN_1になる
#define ST_BRANCH_NM 0
#define EN_BRANCH_NM 50
#define DELTA_BRANCH_NM 10

//屈折率
#define N_0 1.0
#define N_1 1.56
//serikon
//#define N_1 8.4179

//先端における横幅の割合
#define ST_EDGE_RATE 0.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.5

//ラメラの先端を丸める曲率 (0で四角形のまま, 1.0で最もカーブする)
#define CURVE 0.0

static int width_nm[2]     = {ST_WIDTH_NM, ST_WIDTH_NM};
static int thickness_nm[2] = {ST_THICK_NM_0, ST_THICK_NM_1};
static int layerNum = ST_LAYER_NUM;     //枚数
static int branch_width_nm = ST_BRANCH_NM; //枝の幅

static double branch_width_s; //枝の幅
static double width_s[2];     //幅
static double thickness_s[2]; //厚さ

static double ep_s[2];        //誘電率 = n*n*ep0

static double edge_width_rate = ST_EDGE_RATE;

static double c0, c1; //2次関数の比例定数

static double calc_width(double sx, double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
  double new_wid = (wid+branch_width_s)*(p + (1-p)*edge_width_rate);
  //  double new_wid = wid-branch_width_s*2 + (p + (1-p)*edge_width_rate)*branch_width_s;

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
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/split; //細分化したセルの位置
      double sy = _y + row*j/split;

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
      
      if(fabs(sx) < wid/2) //width_s[k]
        s[k] +=1;
    }    
  }

  s[0] /= split*split;
  s[1] /= split*split;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep_s[0]*s[0] + ep_s[1]*s[1];
}

double ( *multiLayerModel_EPS(void))(double, double, int, int)
{
  return eps;
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
  return false;  
}

bool multiLayerModel_isFinish(void)
{
  return nextStructure();
}

void multiLayerModel_needSize(int *x_nm, int *y_nm)
{
  (*x_nm) = max( width_nm[0], width_nm[1]) + branch_width_nm;
  (*y_nm) = (thickness_nm[0]+thickness_nm[1])*layerNum;
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
  sprintf(buf,"n_%.2lf_%.2lf", N_0, N_1);
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

void multiLayerModel_init()
{
  width_s[0]     = field_toCellUnit(width_nm[0]);
  width_s[1]     = field_toCellUnit(width_nm[1]);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[1]);
  ep_s[0] = N_0*N_0*EPSILON_0_S;
  ep_s[1] = N_1*N_1*EPSILON_0_S;

  branch_width_s = field_toCellUnit(branch_width_nm);
  
  c0 = -4*width_s[0]*CURVE/thickness_s[0]/thickness_s[0];
  c1 = -4*width_s[1]*CURVE/thickness_s[1]/thickness_s[1];
}
