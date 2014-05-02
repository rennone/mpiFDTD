#ifndef NTFF_TM_H
#define NTFF_TM_H

#include "myComplex.h"

extern void ntffTM_Frequency( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resEz[360]);
extern void ntffTM_FrequencySplit( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resEz[360]);

extern void ntffTM_TimeCalc(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);
extern void ntffTM_TimeCalc2(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz);
extern void ntffTM_TimeTranslate(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, dcomplex *Eth, dcomplex *Eph);
#endif
