#include "ntffTE.h"
#include "field.h"
#include "function.h"
#include "bool.h"
#include <math.h>
#include <stdlib.h>

static double R0;

void ntffTE_init()
{
  R0 = 1.0e6 * field_toCellUnit(500);//* field_getLambda_S();  
}
//---------------------- ntff--------------------//

void ntffTE_TimeTranslate(dcomplex *Wx, dcomplex *Wy,dcomplex *Uz, dcomplex *Eth, dcomplex *Eph)
{
  const double w_s = field_getOmega();
  // 1/4πRcをここでかけている(宇野先生の本では無限遠方点のためRを省略している)
  // たぶん,Rをつける意味はないと思う. => kesita
  const double complex coef = 1.0/(4*M_PI*C_0_S)*csqrt( 2*M_PI*C_0_S/(I*w_s) );
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();

  double theta = 0;
  //NTFF TE Translate
  for(int ang=0; ang<360; ang++)
  {
    double phi = ang*M_PI/180.0;
    
    int k= ang*nInfo.arraySize;    
    double sx = cos(theta)*cos(phi);
    double sy = cos(theta)*sin(phi);
    double sz = -cos(theta); //宇野先生の本では -sin(theta)になってる
    double px = -sin(phi);
    double py = cos(phi);
    int i;
    for(i=0; i < maxTime; i++)
    {
      double complex WTH = Wx[k+i]*sx + Wy[k+i]*sy + 0;
      double complex WPH = Wx[k+i]*px + Wy[k+i]*py;
      double complex UTH = 0 + 0 + Uz[k+i]*sz;
      double complex UPH = 0 + 0;
      double complex ETH  = coef*(-Z_0_S*WTH-UPH);
      double complex EPH  = coef*(-Z_0_S*WPH+UTH);
      
      Eth[k+i] = ETH;
      Eph[k+i] = EPH;
    }
  }  
}

static inline void calc(double time_plus_timeShift, dcomplex eh,  dcomplex *UW_ang){
  int m = floor(time_plus_timeShift+0.5);
  double a = (0.5 + time_plus_timeShift - m);
  double b = 1.0-a;
  double ab = a-b;
//  const double coef = 1.0/(4*M_PI*C_0_S*R0);
  UW_ang[m-1] += eh*b;//*coef;
  UW_ang[m]   += eh*ab;//*coef;
  UW_ang[m+1] -= eh*a;//*coef;  
}

 void ntffTE_TimeCalc(dcomplex *Ex, dcomplex *Ey, dcomplex *Hz,
                            dcomplex *Wx, dcomplex *Wy,dcomplex *Uz)
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();

  double timeE = field_getTime() - 1;  //t - Δt  todo
  double timeH = field_getTime() - 0.5;  //t - Δt/2  todo

  // 中心から各頂点へのベクトルを計算しておく
  double lt_cx = nInfo.left - nInfo.cx;
  double rt_cx = nInfo.right - nInfo.cx;
  double bm_cy = nInfo.bottom - nInfo.cy;
  double tp_cy = nInfo.top - nInfo.cy;

  //各頂点のインデックス(Ex,Ey,Hz用)も計算しておく
  int lb_index   = field_index(nInfo.left , nInfo.bottom); 
  int lt_index   = field_index(nInfo.left , nInfo.top);
  int rb_index   = field_index(nInfo.right, nInfo.bottom);
  int rt_index   = field_index(nInfo.right, nInfo.top);

  double ToRad = M_PI/180.0;
  for(int ang=0, stp=0; ang<360; ang++, stp+=nInfo.arraySize)
  {
    double rad = ang*ToRad;
    double r1x_per_c = cos(rad)/C_0_S, r1y_per_c = sin(rad)/C_0_S;
    
    //ang°の位置にシフトしたポジション, こうすれば Ux_ang[i]でその角度のi番目にアクセスできる.
    dcomplex *Wx_ang = &Wx[stp];
    dcomplex *Wy_ang = &Wy[stp];
    dcomplex *Uz_ang = &Uz[stp];

    //(l,b)->(r,b) n is (0,-1, 0)  Js = n × H = (-Hz, 0, 0)  Ms = E × n = (0, 0, -Ex)
    {
      const double r2x = lt_cx+0.5, r2y = bm_cy;        
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC;
      for(int k=lb_index; k<rb_index; k+=fInfo_s.N_PY)
      {
        double complex ex = -Ex[k];
        double complex hz = -0.5*( Hz[k] + Hz[k-1] );
        calc(timeE+timeShift, ex, Uz_ang);
        calc(timeH+timeShift, hz, Wx_ang);
        timeShift -= r1x_per_c; //右に1セル進むと -r1x_per_c大きくなる(timeShiftの式より)        
      }
    }

    //(r,b)->(r,t)  n(1,0) Js = n × H = (0,-Hz,0)    Ms = E × n = (0,0,-Ey)
    {
      const double r2x = rt_cx, r2y = bm_cy+0.5;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC;
      for(int k=rb_index; k<rt_index; k++)
      {
        double complex ey = -Ey[k];
        double complex hz = -0.5*( Hz[k] + Hz[k-fInfo_s.N_PY] );
        calc(timeE+timeShift, ey, Uz_ang);
        calc(timeH+timeShift, hz, Wy_ang);
        timeShift -= r1y_per_c; //右に1セル進むと -r1x_per_c大きくなる(timeShiftの式より)        
      }
    }

    //(l,t)->(r,t) n is (0,1)    //Js = n × H = (Hz, 0, 0)  Ms = E × n = (0, 0, Ex)
    {
      const double r2x = lt_cx+0.5, r2y = tp_cy;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC;      
      for ( int k=lt_index; k<rt_index; k+=fInfo_s.N_PY ) 
      {
        double complex ex = Ex[k];
        double complex hz = 0.5*( Hz[k] + Hz[k-1] );

        calc(timeE+timeShift, ex, Uz_ang);
        calc(timeH+timeShift, hz, Wx_ang);
        timeShift -= r1x_per_c; //右に1セル進むと -r1x_per_c大きくなる(timeShiftの式より)        
      }
    }

    //(l,b)->(l,t) n is (-1,0)    //Js = n × H = (0,Hz,0)    Ms = E × n = (0,0,Ey)    
    {
      const double r2x = lt_cx, r2y = bm_cy+0.5;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC;       
      for(int k=lb_index; k<lt_index; k++ )
      {
        double complex ey = Ey[k];
        double complex hz = 0.5*( Hz[k] + Hz[k-fInfo_s.N_PY] );
        calc(timeE+timeShift, ey, Uz_ang);
        calc(timeH+timeShift, hz, Wy_ang);
        timeShift -= r1y_per_c; //右に1セル進むと -r1x_per_c大きくなる(timeShiftの式より)        
      }
    }    
  }
}


void ntffTE_TimeOutput(dcomplex *Wx, dcomplex *Wy, dcomplex *Uz, FILE *fpRe, FILE *fpIm)
{
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();
  dcomplex *Eth, *Eph;  
  Eth = newDComplex(360*nInfo.arraySize);
  Eph = newDComplex(360*nInfo.arraySize);
  ntffTE_TimeTranslate(Wx,Wy,Uz,Eth,Eph);
  for(int ang=0; ang<360; ang++)
  {
    int k = ang*nInfo.arraySize;
    for(int i=0; i < maxTime; i++)
    {
      fprintf(fpRe,"%.20lf " , creal(Eph[k+i]));
      fprintf(fpIm,"%.20lf " , cimag(Eph[k+i]));   
    }
    
    fprintf(fpRe,"\n");
    fprintf(fpIm,"\n");
  }  
  free(Eth);
  free(Eph);
}

//積分路では, セル
//なるべく,補間(平均)につかう値が少ないような位置となる積分路を選ぶ.
void ntffTE_Frequency(dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resultEphi[360])
{
  FieldInfo  fInfo         = field_getFieldInfo();
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();
  WaveInfo_S waveInfo_s    = field_getWaveInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();
  
  double cx = nInfo.cx;
  double cy = nInfo.cy;
  int tp = nInfo.top;
  int bm = nInfo.bottom;
  int rt = nInfo.right;
  int lt = nInfo.left;
  
  dcomplex coef = csqrt( I*waveInfo_s.K_s/(8*M_PI*R0) ) * cexp(I*waveInfo_s.K_s*R0);  

  //Js -> W   Ms -> U
  for ( int ang=0; ang<360; ang++ )
  {
    double rad = ang * M_PI/180.0;
    double rx = cos(rad), ry = sin(rad);
    
    dcomplex Nx = 0;
    dcomplex Ny = 0;
    dcomplex Lz = 0;

    //(l,b)->(r,b) n(0,-1,0)  Js = n × H = (-Hz, 0, 0)  Ms = E × n = (0, 0, -Ex)
    for( int i=lt; i<=rt; i++)
    {
      const double r2x = i-cx+0.5;
      const double r2y = bm-cy;

      int k = field_index(i, bm);
      double complex C_HZ  = 0.5*( Hz[k] + Hz[k-1] );
      double complex C_EX  = Ex[k];

      double innerProd = rx*r2x + ry*r2y;
      Nx   -= C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
      Lz   -= C_EX * cexp( I * waveInfo_s.K_s * innerProd );
    }

    //(r,b)->(r,t)  n(1,0) Js = n × H = (0,-Hz,0)    Ms = E × n = (0,0,-Ey)
    for(int j=bm; j<tp; j++)
    {
      double r2x = rt-cx;
      double r2y = j-cy+0.5;  

      int k = field_subIndex(rt, j);
      double complex C_HZ  = 0.5*( Hz[k] + Hz[k-fInfo_s.N_PY] );
      double complex C_EY  = Ey[k];

      double innerProd = rx*r2x + ry*r2y;  //内積
      Ny -= C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
      Lz -= C_EY * cexp( I * waveInfo_s.K_s * innerProd );
    }


    //(l,t)->(r,t) n is (0,1)    //Js = n × H = (Hz, 0, 0)  Ms = E × n = (0, 0, Ex)
    for ( int i=lt; i<rt; i++ ) 
    {
      const double r2x = i-cx+0.5;
      const double r2y = tp-cy;

      int k = field_subIndex(i, tp);
      double complex C_HZ  = 0.5*( Hz[k]+ Hz[k-1] );
      double complex C_EX  = Ex[k];
      double innerProd = rx*r2x  + ry*r2y;  //内積
      Nx += C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
      Lz += C_EX * cexp( I * waveInfo_s.K_s * innerProd );
    }
    
    //(l,b)->(l,t) n is (-1,0)    //Js = n × H = (0,Hz,0)    Ms = E × n = (0,0,Ey)    
    for ( int j=bm; j<tp; j++ )
    {
      const double r2x = lt-cx;
      const double r2y = j-cy+0.5;
        
      int k = field_subIndex(lt, j);
      double complex C_HZ  = 0.5 * ( Hz[k] + Hz[k-fInfo_s.N_PY] );
      double complex C_EY  = Ey[k];
      double innerProd = rx*r2x  + ry*r2y;  //内積
      Ny += C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
      Lz += C_EY * cexp( I * waveInfo_s.K_s * innerProd );
    }    
    
    // Get Ephi
    double complex Nphi  = -Nx*sin(rad) + Ny*cos(rad);
    resultEphi[ang] = coef * ( Z_0_S*Nphi - Lz )*fInfo.h_u_nm; //物理単位に変換しておく
  }

  FILE *fp = fopen("ntffStrTE.txt", "w");  
  for(int ang = 0; ang<360; ang++)
  {
    fprintf(fp, "%.18lf\n", cnorm(resultEphi[ang]));
  }
  printf("ntff te frequency finished\n");
}

//--------------------------------------分割領域バージョン------------------------------------------
void ntffTE_FrequencySplit(dcomplex *Ex, dcomplex *Ey, dcomplex *Hz, dcomplex resultEphi[360])
{
  FieldInfo  fInfo         = field_getFieldInfo();
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  WaveInfo_S waveInfo_s    = field_getWaveInfo_S();
  
  //分割領域系に置ける, 中心の位置
  double cx = fInfo_s.N_PX/2 - subInfo_s.OFFSET_X;
  double cy = fInfo_s.N_PY/2 - subInfo_s.OFFSET_Y;

  dcomplex coef = csqrt( I*waveInfo_s.K_s/(8*M_PI*R0) ) * cexp(I*waveInfo_s.K_s*R0);
  
  NTFFInfo nInfo = field_getNTFFInfo();
  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    nInfo.top - subInfo_s.OFFSET_Y; //上面
  int bm = nInfo.bottom - subInfo_s.OFFSET_Y; //下面
  int rt = nInfo.right  - subInfo_s.OFFSET_X; //右
  int lt = nInfo.left   - subInfo_s.OFFSET_X; //左

  //Js -> W                 Ms -> U
  for ( int ang=0; ang<360; ang++ )
  {
    double rad = ang * M_PI/180.0;
    double rx = cos(rad), ry = sin(rad);
    
    dcomplex Nx = 0;
    dcomplex Ny = 0;
    dcomplex Lz = 0;
    
    //bottom side
    //normal vector n is (0,-1, 0)   //Js = n × H = (-Hz, 0, 0)  Ms = E × n = (0, 0, -Ex)
    if ( 0 < bm && bm < subInfo_s.SUB_N_PY-1 )
    {
      //積分路の端で無ければ, 分割領域の端から端までが積分路である
      const int subLeft  = max(1, lt);
      const int subRight = min(subInfo_s.SUB_N_PX-1, rt);
      int i;
      for( i=subLeft; i<subRight; i++)
      {
        const double r2x = i-cx+0.5;
        const double r2y = bm-cy;            //distance between origin to cell

        int k = field_subIndex(i, bm);
        double complex C_HZ  = 0.5*( Hz[k] + Hz[k-1] );
        double complex C_EX  = Ex[k];

        double innerProd = rx*r2x + ry*r2y;
        Nx   -= C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
        Lz   -= C_EX * cexp( I * waveInfo_s.K_s * innerProd );
      }
    }
      //right side
    //normal vector n is (1,0)
    //Js = n × H = (0,-Hz,0)    Ms = E × n = (0,0,-Ey)
    if ( 0 < rt && rt < subInfo_s.SUB_N_PX-1)
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm);
      int j;
      for(j=subBottom; j<subTop; j++)
      {
        double r2x = rt-cx;
        double r2y = j-cy+0.5;  

        int k = field_subIndex(rt, j);
        double complex C_HZ  = 0.5*( Hz[k] + Hz[k-subInfo_s.SUB_N_PY] );
        double complex C_EY  = Ey[k];

        double innerProd = rx*r2x + ry*r2y;  //内積
        Ny -= C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
        Lz -= C_EY * cexp( I * waveInfo_s.K_s * innerProd );
      }
    }    


    //top side
    //normal vector n is (0,1)    //Js = n × H = (Hz, 0, 0)  Ms = E × n = (0, 0, Ex)
    if ( 0 < tp && tp < subInfo_s.SUB_N_PY-1)
    {      
      int subLeft  = max(1, lt);
      int subRight = min(subInfo_s.SUB_N_PX-1, rt);
      int i;
      for ( i=subLeft; i<subRight; i++ ) 
      {
        const double r2x = i-cx+0.5;
        const double r2y = tp-cy;

        int k = field_subIndex(i, tp);
        double complex C_HZ  = 0.5*( Hz[k]+ Hz[k-1] );
        double complex C_EX  = Ex[k];
        double innerProd = rx*r2x  + ry*r2y;  //内積
        Nx += C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
        Lz += C_EX * cexp( I * waveInfo_s.K_s * innerProd );
      }
    }

    //left side
    //normal vector n is (-1,0)    //Js = n × H = (0,Hz,0)    Ms = E × n = (0,0,Ey)    
    if ( 0 < lt && lt < subInfo_s.SUB_N_PX-1 )
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm);
      int j;
      for ( j=subBottom; j<subTop; j++ )
      {
        const double r2x = lt-cx;
        const double r2y = j-cy+0.5;
        
        int k = field_subIndex(lt, j);
        double complex C_HZ  = 0.5 * ( Hz[k] + Hz[k-subInfo_s.SUB_N_PY] );
        double complex C_EY  = Ey[k];
        double innerProd = rx*r2x  + ry*r2y;  //内積
        Ny += C_HZ * cexp( I * waveInfo_s.K_s * innerProd );
        Lz += C_EY * cexp( I * waveInfo_s.K_s * innerProd );
      }
    }
    
    // Get Ephi
    double complex Nphi  = -Nx*sin(rad) + Ny*cos(rad);
    resultEphi[ang] = coef * ( Z_0_S*Nphi - Lz )*sqrt(fInfo.h_u_nm);//物理単位に変換する
  }

  FILE *fp = fopen("ntffStrTE.txt", "w");  
  for(int ang = 0; ang<360; ang++)
  {
    fprintf(fp, "%.18lf\n", cnorm(resultEphi[ang]));
  }
  printf("ntff te frequency finished\n");
}
