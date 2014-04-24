#include "ntffTM.h"
#include "field.h"
#include <math.h>
/*
void ntffTM_Frequency( double complex *Hx, double complex *Hy, double complex *Ez, double complex resEz[360])
{
  double r0 = 1.0e6;
    
  const double cx = N_PX/2 - offsetX;  //分割領域系に置ける, 中心の位置
  const double cy = N_PY/2 - offsetY;
  
  double k_s = field_getK();
  double complex coef = csqrt( I*k_s/(8*M_PI*r0) ) * cexp(I*k_s*r0);  
  NTFFInfo nInfo = field_getNTFFInfo();
  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    nInfo.top - offsetY; //上面
  int bm = nInfo.bottom - offsetY; //下面
  int rt = nInfo.right  - offsetX; //右
  int lt = nInfo.left   - offsetX; //左

  const int max_angle = 360;
  int ang;
  static double complex ntffEz[360];

  for(ang=0; ang<max_angle; ang++)
  {
    double rad = ang*M_PI/180.0;
    double rx  = cos(rad), ry = sin(rad);
    double r2x, r2y;

    double complex Nz = 0;
    double complex Lx = 0;
    double complex Ly = 0;
    double complex C_EZ, C_HX, C_HY;
		
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    //積分路が分割領域内にあるか確認
    if (  0 < bm && bm < SUB_N_PY-1)
    {
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで
      for ( int i=subLeft; i<subRight; i++ )
      {
        r2x  =  i-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        r2y  = bm-cy;

        int k = subInd(&i, &bm);
        C_EZ = Ez[k];
        C_HX = 0.5 * ( Hx[k] + Hx[subIndBottom(k)] );
        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0   ,0,C_HX)
        Lx  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (C_EZ,0,   0)
      }
    }
    
    // (right,bottom) -> (right,top) n=(1,0)
    if ( 0 < rt && rt < SUB_N_PX-1)
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ )
      {
        r2x  = rt-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        r2y  =  j-cy;

        int k = subInd(&rt, &j);
        C_EZ = Ez[k];
        C_HY = 0.5 * ( Hy[k] + Hy[subIndLeft(k)] );

        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,   0, C_HY)
        Ly  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0,C_EZ,    0)
      }
    }

    // (right,top) -> (left,top)  n=(0,1)
    if ( 0 < tp && tp < SUB_N_PY-1)
    {      
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで      
      for ( int i=subLeft; i<subRight; i++ )
      {
        r2x  =  i-cx;
        r2y  = tp-cy;

        int k = subInd(&i, &tp);
        C_EZ = Ez[k];
        C_HX = 0.5 * ( Hx[k] + Hx[subIndBottom(k)] );

        double innerProd = rx*r2x  + ry*r2y;  //内積
        Nz   -= C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0,    0, -C_HX)
        Lx   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
      }
    }

    // (left,top) -> (left,bottom)
    if ( 0 < lt && lt < SUB_N_PX )
    {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm); 
      for ( int j=subBottom; j<subTop; j++ )
      {
        r2x =   lt-cx; r2y = j-cy;

        int k = subInd(&lt, &j);
        C_EZ  = Ez[k];
        C_HY  = 0.5 * ( Hy[k] + Hy[subIndLeft(k)] );

        double innerProd = rx*r2x  + ry*r2y;  //内積      
        Nz   -= C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
        Ly   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
      }
    }

    double complex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    ntffEz[ang] = coef * ( Z_0_S*Nz + Lphi );
  }

  const char nameRe[256]="ntffRe.txt";
  const char nameIm[256]="ntffIm.txt";
  FILE *fpR = fileOpen(nameRe);
  FILE *fpI = fileOpen(nameIm);
  
  for(ang = 0; ang<max_angle; ang++)
  {
    fprintf(fpR, "%.18lf\n", creal(ntffEz[ang]));
    fprintf(fpI, "%.18lf\n", cimag(ntffEz[ang]));
  }

  fclose(fpR);
  fclose(fpI);
  printf("saved at %s & %s \n", nameRe, nameIm);
  printf("ntffFrequency end\n");
}
*/
