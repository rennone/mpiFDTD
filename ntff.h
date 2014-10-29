#ifndef __NTFF_H
#define __NTFF_H

#define LAMBDA_ST_NM 380
#define LAMBDA_EN_NM 700
#define NTFF_NUM 8192 // 1<<13

extern void ntff_outputEnormTxt(double **e_norm, const char *file_name);
extern void ntff_outputEnormBin(double **e_norm, const char *file_name);

#endif
