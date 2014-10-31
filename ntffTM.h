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

//時間領域の遠方解 Ux, Uy, WzからEth, Ephを求め, fpRe, fpImに書き出す
extern void ntffTM_TimeOutput(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);

#endif
