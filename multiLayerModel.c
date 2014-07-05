#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラの厚さ
#define ST_THICK_NM 90
#define EN_THICK_NM 160
#define DELTA_THICK_NM 10

//ラメラの枚数
#define LAYER_NUM 8

//互い違い
#define ASYMMETRY false

//屈折率
#define N_1 1.0
#define N_2 1.56

static int width_nm[2] = {ST_WIDTH_NM, ST_WIDTH_NM};
static int thickness_nm[2] = {ST_THICK_NM, ST_THICK_NM};
static int layerNum = LAYER_NUM;     //枚数

static double width_s[2];     //幅
static double thickness_s[2]; //厚さ
static double ep_s[2];        //誘電率 = n*n*ep0

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
      
      if(abs(sx) < width_s[k]/2)
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
  return thickness_nm[0] > EN_THICK_NM;
}

void multiLayerModel_needSize(int *x_nm, int *y_nm)
{
  (*x_nm) = max( width_nm[0], width_nm[1]);
  (*y_nm) = (thickness_nm[0]+thickness_nm[1])*LAYER_NUM;
}

void multiLayerModel_moveDirectory()
{
  char buf[512];
  sprintf(buf, "thick_%dnm_layer_%d", thickness_nm[0], layerNum);
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
}
