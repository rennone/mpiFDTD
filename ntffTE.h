#ifndef NTFF_TE_H
#define NTFF_TE_H
#include "myComplex.h"
#include <stdio.h>
extern void ntffTE_Frequency( dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resEphi[360]);
extern void ntffTE_FrequencySplit( dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resEphi[360]);

extern void ntffTE_TimeCalc(dcomplex *Ex, dcomplex *Ey, dcomplex *Hz,
                            dcomplex *Wx, dcomplex *Wy, dcomplex *Uz);

extern void ntffTE_TimeTranslate(dcomplex *Wx, dcomplex *Wy, dcomplex *Uz,dcomplex *Eth, dcomplex *Eph);

//時間領域の遠方解 Wx, Wy, UzからEth, Ephを求め, fpRe, fpImに書き出す
extern void ntffTE_TimeOutput(dcomplex *Wx, dcomplex *Wy, dcomplex *Uz, FILE *fpRe, FILE *fpIm);

#endif
