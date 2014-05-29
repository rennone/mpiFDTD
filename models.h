#ifndef _MODELS_H
#define _MODELS_H

enum MODEL
{
  NO_MODEL,
  MIE_CYLINDER,
  SHELF,
  NONSHELF,
  LAYER,
  MORPHO_SCALE
};

enum MODE{
  D_X, //x方向に線積分
  D_Y, //y方向に線積分
  D_XY //面積分 
};

//モデルを定義
extern void setModel(enum MODEL model);

//(x,y)における誘電率を取得
extern double models_eps(double x, double y, enum MODE mode);
extern bool mpdels_isFinish(void);
extern void models_moveDirectory(void); //所定のディレクトリまで移動する為の関数
//データの吐き出し
//extern void models_output(FILE *fp, double complex *data);
#endif
