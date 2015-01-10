#ifndef NTFF_TM_H
#define NTFF_TM_H

#include <stdio.h>
#include "myComplex.h"

extern void ntffTM_init();

//resEzに周波数領域の遠方解を代入する.
extern void ntffTM_Frequency( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resEz[360]);

//fpRe, fpImに周波数領域の遠方解を書き出す
extern void ntffTM_FreqOutput(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, FILE *fpRe, FILE *fpIm);

//resEzに周波数領域の遠方解を代入する(領域分割タイプ).
extern void ntffTM_FrequencySplit( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resEz[360]);

//時間領域の遠方解のupdate処理(Ux,Uy,Wzを更新)
extern void ntffTM_TimeCalc(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);

//時間領域の遠方解のupdate処理(Ux,Uy,Wzを更新)
//extern void ntffTM_TimeCalc2(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);

//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求める.
extern void ntffTM_TimeTranslate(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, dcomplex *Eth, dcomplex *Eph);

//評価(遺伝的アルゴリズム)に反射率を用いる為, 反射率の配列を返したい
//なので, 下にあるntffTM_TimeCalcNormを用いたほうが良い
//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求め, fpRe, fpImに書き出す
extern void ntffTM_TimeOutput(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);

//時間領域の遠方解 Ux, Uy, WzからEthを求め
//それをfftした結果から波長,反射角ごとのEthのノルムを返す
// res[λ][deg]でアクセスできる. 0 <= deg <= 359, stLambda <= λ <= enLamba
extern double** ntffTM_TimeCalcNorm(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, int *stLambda, int *enLambda);
#endif
