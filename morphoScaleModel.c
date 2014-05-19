#include "morphoScaleModel.h"
#include "field.h"
#include "function.h"
#include <math.h>

double width_s[2];     //幅
double widthOnTopRate; //頂上だと幅が何倍になるか(基本0~1)
double thickness_s[2]; //厚さ
double ep[2];           //誘電率 = n*n*ep0
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
  return EPSILON_0_S*(1-s[0]-s[1]) + ep[0]*s[0] + ep[1]*s[1];
}

#include <stdlib.h>

double ( *morphoScaleModel_EPS(void))(double, double, int, int)
{
  FILE *fp = NULL;
  if( !(fp = fopen("config.txt", "r")))
  {
    printf("cannot find config.txt of morphoScaleModel\n");
    exit(2);
  }
  int err;
  char buf[1024];
  while( fgets(buf, 1024, fp) != NULL)
  {
    if(buf[0] == '#' || buf[0] == '\0' || buf[0] == '\n')
      continue;
    // #より後ろはカット
    char* p = strstr(buf,"#"); //#の位置を探す
    if( p != NULL){
      strncpy(buf, buf, p-buf+1);
    }
    printf("%d\n",atoi(buf));
    //fscanf(fp, "width=%lf;%lf")
  }
  width_s[0]     = field_toCellUnit(300);
  width_s[1]     = field_toCellUnit(300);
  thickness_s[0] = field_toCellUnit(90);
  thickness_s[1] = field_toCellUnit(90);
  layerNum = 8;
  ep[0] = 1.56*1.56*EPSILON_0_S;
  ep[1] = 1*1*EPSILON_0_S;

  widthOnTopRate = 0.2; //頂上になると幅が半分になる
  notsymmetry = false;
  return eps;
}
