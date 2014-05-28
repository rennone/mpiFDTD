#include "multiLayerModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

double width_s[2];     //幅
double thickness_s[2]; //厚さ
double ep[2];           //誘電率 = n*n*ep0
int layerNum;          //枚数
bool notsymmetry;      //左右比対称

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
      
      if(abs(sx) < width_s[k]/2)
        s[k] +=1;     

    }    
  }

  s[0] /= 32.0*32.0;
  s[1] /= 32.0*32.0;
  return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

double ( *multiLayerModel_EPS(void))(double, double, int, int)
{
  width_s[0]     = field_toCellUnit(300);
  width_s[1]     = field_toCellUnit(300);
  thickness_s[0] = field_toCellUnit(90);
  thickness_s[1] = field_toCellUnit(90);
  layerNum = 8;
  ep[0] = 1.56*1.56*EPSILON_0_S;
  ep[1] = 1*1*EPSILON_0_S;

  notsymmetry = true;

  return eps;
}

bool multiLayerModel_isFinish(void)
{
  return true;
}
