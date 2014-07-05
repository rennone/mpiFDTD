#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラの厚さ
#define ST_THICK_NM 30
#define EN_THICK_NM 160
#define DELTA_THICK_NM 10

//ラメラの枚数
#define LAYER_NUM 8

//互い違い
#define ASYMMETRY true

//屈折率
#define N_1 1.0
#define N_2 1.56

static int width_nm[2] = {ST_WIDTH_NM, ST_WIDTH_NM};
static int thickness_nm[2] = {ST_THICK_NM, ST_THICK_NM};
static int layerNum = LAYER_NUM;     //枚数

static double width_s[2];     //幅
static double thickness_s[2]; //厚さ

static double ep_s[2];        //誘電率 = n*n*ep0

//先端における横幅の割合

#define ST_EDGE_RATE 1.0
#define EN_EDGE_RATE 1.0
#define DELTA_EDGE_RATE 0.1
static double edge_width_rate = ST_EDGE_RATE;

//ラメラの曲率 (1で四角形のまま, 0.0で最もカーブする)
#define CURVE 1.0
static double c;

static double calc_width(double sy, double wid, double hei, double modY, int k)
{
  double p = 1 - sy/hei;
  double new_wid = wid*(p + (1-p)*edge_width_rate);

  //2次関数で曲率計算
  double dh = k==0 ? modY : modY - thickness_s[0];
  
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
  double split = 32;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = _x + col*i/split; //細分化したセルの位置
      double sy = _y + row*j/split;

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

      if (sx < 0 && ASYMMETRY)
        k = 1-k;		//左右で反転, 互い違いでなかったら反転しない

//      double p = 1 - sy/height;
      double wid = calc_width(sy, width_s[k], height, modY, k);
      if(abs(sx) < wid/2) //width_s[k]
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

bool multiLayerModel_isFinish(void)
{
  thickness_nm[0] += DELTA_THICK_NM;
  thickness_nm[1] += DELTA_THICK_NM;

  if(thickness_nm[0] > EN_THICK_NM)
  {
    thickness_nm[0] = ST_THICK_NM;
    thickness_nm[1] = ST_THICK_NM;

    edge_width_rate += DELTA_EDGE_RATE;
    if(edge_width_rate > EN_EDGE_RATE)
      return true;
  }
     
  return false;
}

void multiLayerModel_needSize(int *x_nm, int *y_nm)
{
  (*x_nm) = max( width_nm[0], width_nm[1]);
  (*y_nm) = (thickness_nm[0]+thickness_nm[1])*LAYER_NUM;
}

void multiLayerModel_moveDirectory()
{
  char buf[512];
  sprintf(buf, "thick_%dnm_layer_%d_sym_%d", thickness_nm[0], layerNum, ASYMMETRY);
  makeDirectory(buf);
  moveDirectory(buf);
}

void multiLayerModel_init()
{
  width_s[0]     = field_toCellUnit(width_nm[0]);
  width_s[1]     = field_toCellUnit(width_nm[1]);
  thickness_s[0] = field_toCellUnit(thickness_nm[0]);
  thickness_s[1] = field_toCellUnit(thickness_nm[1]);
  ep_s[0] = N_1*N_1*EPSILON_0_S;
  ep_s[1] = N_2*N_2*EPSILON_0_S;

  c = 4*width_s[0]*(CURVE-1)/thickness_s[0]/thickness_s[0];
}
