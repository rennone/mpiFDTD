#include "zigzagModel.h"
#include "function.h"
#include "field.h"
#include <math.h>

//横幅
#define ST_WIDTH_NM 300
#define EN_WIDTH_NM 300
#define DELTA_WIDTH_NM 10

//ラメラの厚さ
#define ST_THICK_NM 30
#define EN_THICK_NM 150
#define DELTA_THICK_NM 30

//ラメラの枚数
#define ST_LAYER_NUM 11
#define EN_LAYER_NUM 11
#define DELTA_LAYER_NUM 1

//屈折率
#define N_0 1.56

#define ST_DEGREE 10
#define EN_DEGREE 80
#define DELTA_DEGREE 10

static int width_nm     = ST_WIDTH_NM;
static int thickness_nm = ST_THICK_NM;
static int layerNum     = ST_LAYER_NUM;     //枚数

static int degree = ST_DEGREE;
static double rad;
static double width_s;     //幅
static double thickness_s; //厚さ
static double height_s;    //全体の高さ.
static double ox_s, oy_s;  //原点
static double expand_width_s, expand_height_s; //はみ出る分

static double ep_s;        //誘電率 = n*n*ep0

static double eps(double _x, double _y, int col, int row)
{
  double x = _x - ox_s;
  double y = _y - oy_s;
 
  if( fabs(x) > width_s/2.0 + expand_width_s + 0.5 || y < -expand_height_s-0.5 || y > height_s+expand_height_s+0.5)
    return EPSILON_0_S;

  double s=0; //n1の分割セルの数が入る  
  double split = 10;
  double half_split = split/2;
  for(double i=-half_split+0.5; i<half_split; i+=1){
    for(double j=-half_split+0.5; j<half_split; j+=1){
      double sx = x + col*i/split; //細分化したセルの位置
      double sy = y + row*j/split;

      if( fabs(sx) > width_s/2 + expand_width_s || sy < -expand_height_s || sy > height_s+expand_height_s)
        continue;

      //deg = 0の時に,線が縦に積まれるように(重ならないように),cos*thick/2だけずらしている.
      int k = min( layerNum-1, max( 0, sy/(width_s*sin(rad) + cos(rad)*thickness_s/2) ));
      //直線の式
      double px = -width_s*cos(rad)/2.0;          //位置ベクトル
      double py = floor((k+1)/2)*2*(width_s*sin(rad)) + k*cos(rad)*thickness_s/2;
      double vx = (1 - ((k&1)<<1))*cos(rad); //方向ベクトル,偶奇で逆になる
      double vy = sin(rad);

      //直線との距離を求める
      double sqrDist = pow(sx-px,2)+pow(sy-py,2) - pow((sx-px)*vx + (sy-py)*vy, 2);

      if(sqrDist <= thickness_s*thickness_s)
        s+=1;
    }
  }  
  s /= split*split;
  return EPSILON_0_S*(1-s) + ep_s*s;
}

static bool nextStructure()
{
  UN_DONE("zigzag nextStructure");
  return true;
}

double ( *zigzagModel_EPS())(double, double, int, int)
{
  return eps;
}

bool zigzagModel_isFinish()
{
  return nextStructure();
}

void zigzagModel_moveDirectory()
{
  makeDirectory("tmp");
  moveDirectory("tmp");
}

void zigzagModel_needSize(int *x, int*y)
{
  //角度が違うと大きくサイズが異なり, image画像が角度が違うだけに見えないので
  //パフォーマンスは落ちるけど, あえて常に最大サイズの計算領域を確保させる.
  double max_rad = EN_DEGREE * M_PI / 180.0;
  *x = cos(max_rad)*width_nm + thickness_nm*sin(max_rad);
  *y = sin(max_rad)*width_nm*layerNum + (1+layerNum)*thickness_nm*cos(max_rad);
}

void zigzagModel_init()
{
  rad = degree * M_PI / 180.0;
  width_s         = field_toCellUnit(width_nm);
  thickness_s     = field_toCellUnit(thickness_nm);
  height_s        = (width_s*sin(rad))*layerNum + layerNum*thickness_s/2*cos(rad);
  expand_width_s  = thickness_s/2.0*sin(rad);
  expand_height_s = thickness_s/2.0*cos(rad);

  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  ox_s = fInfo_s.N_PX/2.0;
  oy_s = (fInfo_s.N_PY - height_s)/2.0;

  ep_s = N_0*N_0*EPSILON_0_S;
}
