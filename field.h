#ifndef _FIELD_H
#define _FIELD_H
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "bool.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

//入射波のモード
enum WAVE_MODE{
  POINT_LIGHT_IN_CENTER,  //中心に点光源
  SCATTER //散乱波
};

typedef struct NFFInfo
{
  int top, bottom, left, right;
  int cx,cy;
  double RFperC; 
  int arraySize; //必要な配列サイズ
} NTFFInfo;

//シミュレーション上の物理定数
#define C_0_S 0.7071 //下の変数名長いからこっちにする
static const double LIGHT_SPEED_S = 0.7071;

static const double EPSILON_0_S = 1.0;
static const double MU_0_S = 1.0/C_0_S/C_0_S;
static const double Z_0_S  = 1.41422712488; // = √(μ/ε);

extern int N_X;
extern int N_Y;
extern int N_CELL;
extern int N_PML;
extern int N_PX;
extern int N_PY;

extern int field_getOffsetX();
extern int field_getOffsetY();
extern int field_getSubNx();
extern int field_getSubNy();
extern int field_getSubNpx();
extern int field_getSubNpy();
extern int field_getSubNcell();

//インデックスを取ってくる 
extern int ind(const int, const int);

//フィールドの横,縦の大きさ, 1セルのサイズ, pmlレイヤの数, 波長(nm), 計算ステップ
extern void setField(const int wid, const int hei, const double h, const int pml, const double lambda, const double step);

//pml用のσを取ってくる
extern double field_sigmaX(double x, double y);
extern double field_sigmaY(double x, double y);
extern double field_pmlCoef(double x, double y);
extern double field_pmlCoef_LXY(double x, double y);
extern double field_toCellUnit(const double);
extern double field_toPhisycalUnit(const double);

//---------------入射波---------------
extern double complex field_pointLight(void);

//:NTFF
extern void field_nextStep(void);
extern bool field_isFinish(void);


//:getter
extern double field_getT(void);
extern double field_getK(void);
extern double field_getRayCoef(void);
extern double field_getOmega(void);
extern double field_getLambda(void);
extern double field_getWaveAngle(void);
extern double field_getTime(void);
extern double field_getMaxTime(void);
extern NTFFInfo field_getNTFFInfo(void);

//output method
extern void field_outputElliptic(const char *fileName,double complex* data); //
extern void field_outputAllDataComplex(const char *fileName,double complex* data); //
extern void field_outputAllDataDouble(const char *fileName,double* data); //
//--------------------for debug--------------------//
static inline void field_debugPrint(double complex *A)
{
  for(int i=1; i<N_PX; i++){
    for(int j=1; j<N_PY; j++){
      double complex a = A[ind(i,j)];
      if(creal(a) > 1.0)
	printf("(%d, %d) : %lf \n", i, j, creal(a));
    }
  }
}
#endif
