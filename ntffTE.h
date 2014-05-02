#ifndef NTFF_TE_H
#define NTFF_TE_H
#include "myComplex.h"
extern void ntffTE_Frequency( dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resEphi[360]);
extern void ntffTE_FrequencySplit( dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resEphi[360]);

extern void ntffTE_TimeCalc(dcomplex *Ex, dcomplex *Ey, dcomplex *Hz,
                            dcomplex *Wx, dcomplex *Wy, dcomplex *Uz);

extern void ntffTE_TimeTranslate(dcomplex *Wx, dcomplex *Wy, dcomplex *Uz,dcomplex *Eth, dcomplex *Eph);
#endif
