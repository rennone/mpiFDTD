#include <unistd.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "ntff.h"
#include "ntffTM.h"
#include "field.h"
#include "function.h"
#include "bool.h"
#include "cfft.h"

/*
  bottom(k) = k-1
  top(k)    = k+1
  left(k)   = k-SUB_N_PY
  right(k)  = k+SUB_N_PY
*/

static double R0;

static int sub_tp, sub_bm, sub_rt, sub_lt;
static bool IN_TP, IN_BM, IN_LT, IN_RT;
static int sub_ylt, sub_yrt;
static int sub_xtp, sub_xbm;
//static dcomplex *Ux, *Uy, *Wz;

void ntffTM_init()
{
  NTFFInfo nInfo = field_getNTFFInfo();
  
  R0 = 1.0e6 * field_toCellUnit(500);
  
  int tp = nInfo.top;    int bm = nInfo.bottom;  //上下
  int rt = nInfo.right;  int lt = nInfo.left;	 //左右

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  
  sub_tp = tp - subInfo_s.OFFSET_Y;  sub_bm = bm - subInfo_s.OFFSET_Y;
  sub_rt = rt - subInfo_s.OFFSET_X;  sub_lt = lt - subInfo_s.OFFSET_X;

  //以下どれかでも満たせば積分路上に無い
  bool outX = sub_rt <= 0 || sub_lt >= subInfo_s.SUB_N_PX-1; //rtより右, もしくはltより左の小領域
  bool outY = sub_tp <= 0 || sub_bm >= subInfo_s.SUB_N_PY-1; //tpより上, もしくはbmより下の小領域
  
  // 小領域内にどの積分面が存在するか
  IN_TP = (0 < sub_tp && sub_tp < subInfo_s.SUB_N_PY-1) && !outX;
  IN_BM = (0 < sub_bm && sub_bm < subInfo_s.SUB_N_PY-1) && !outX;
  IN_RT = (0 < sub_rt && sub_rt < subInfo_s.SUB_N_PX-1) && !outY;
  IN_LT = (0 < sub_lt && sub_lt < subInfo_s.SUB_N_PX-1) && !outY;

  sub_ylt=-1, sub_yrt=-2;
  if(IN_TP || IN_BM)
  {
    sub_yrt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_rt) );
    sub_ylt = min(subInfo_s.SUB_N_PX-2, max( 1, sub_lt) );
  }

  sub_xtp=-2, sub_xbm=-1;
  if(IN_RT || IN_LT)
  {
    sub_xbm = min(subInfo_s.SUB_N_PY-2, max( 1, sub_bm+1) );  //bm,tpですでに計算しているため, ひとつずれる
    sub_xtp = min(subInfo_s.SUB_N_PY-2, max( 1, sub_tp-1) );  //
  }
}

void ntffTM_finish()
{
  
}

//周波数領域のNTFF
void ntffTM_Frequency( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resultEz[360])
{
  FieldInfo fInfo = field_getFieldInfo();
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();
  double cx = nInfo.cx;
  double cy = nInfo.cy;

  double k_s = field_getK();

  double complex coef = csqrt( I*k_s/(8*M_PI*R0) ) * cexp(I*k_s*R0);
  int tp = nInfo.top;     //上面
  int bm = nInfo.bottom; //下面
  int rt = nInfo.right;  //右
  int lt = nInfo.left;	 //左

  for(int ang=0; ang<360; ang++) {
    double rad = ang*M_PI/180.0;

    double rx  = cos(rad), ry = sin(rad);
    double r2x, r2y;

    dcomplex Nz = 0;
    dcomplex Lx = 0;
    dcomplex Ly = 0;
    dcomplex C_EZ, C_HX, C_HY;
    int i,j;
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    for ( i=lt; i<rt; i++ )
    {
      r2x  =  i-cx;
      r2y  = bm-cy;
      int k = field_index(i,bm);
      C_EZ = Ez[k];
      C_HX = 0.5 * ( Hx[k] + Hx[k-1] );
      
      double innerProd = rx*r2x + ry*r2y;  //内積
      Nz  += C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0   ,0,C_HX)
      Lx  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (C_EZ,0,   0)
    }

    // (right,bottom) -> (right,top) n=(1,0)
    for ( j=bm; j<tp; j++ )
    {
      r2x  = rt-cx;
      r2y  =  j-cy;
      int k = field_index(rt,j);
      C_EZ = Ez[k];
      C_HY = 0.5 * ( Hy[k] + Hy[k-fInfo_s.N_PY] );

      double innerProd = rx*r2x + ry*r2y;  //内積
      Nz  += C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,   0, C_HY)
      Ly  += C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0,C_EZ,    0)
    }

    // (right,top) -> (left,top)  n=(0,1)
    for ( i=lt; i<rt; i++ )
    {
      r2x  =  i-cx;
      r2y  = tp-cy;
      int k = field_index(i,tp);
      C_EZ = Ez[k];
      C_HX = 0.5 * ( Hx[k] + Hx[k-1] );

      double innerProd = rx*r2x  + ry*r2y;  //内積
      Nz   -= C_HX * cexp( I * k_s * innerProd );  // J = n × H = (0,    0, -C_HX)
      Lx   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
    }

    // (left,top) -> (left,bottom)
    for ( j=bm; j<tp; j++ )
    {
      r2x =   lt-cx; r2y = j-cy;
      int k = field_index(lt,j);
      C_EZ  = Ez[k];
      C_HY  = 0.5 * ( Hy[k] + Hy[k-fInfo_s.N_PY]);

      double innerProd = rx*r2x  + ry*r2y;  //内積      
      Nz   -= C_HY * cexp( I * k_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
      Ly   -= C_EZ * cexp( I * k_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
    }

    double complex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    resultEz[ang] = coef * ( Z_0_S*Nz + Lphi )*sqrt(fInfo.h_u_nm);
  }
}


void ntffTM_TimeTranslate(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz, dcomplex *Eth, dcomplex *Eph)
{
  const double w_s = field_getOmega();
  // don't divide by R0 -> textbook by uno
  const double complex coef = 1.0/(4*M_PI*C_0_S)*csqrt( 2*M_PI*C_0_S/(I*w_s) );
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();  
  double theta = 0;
  double ToRad = M_PI/180.0;
  //TM Time Translate
  for(int ang=0, k=0; ang<360; ang++, k+=nInfo.arraySize)
  {
    double phi = ang*ToRad;
    double sx = cos(theta)*cos(phi);
    double sy = cos(theta)*sin(phi);
    double sz = -cos(theta); //宇野先生の本では -sin(theta)になってる(式としては本の方が正しいけどこのプログラムではこうしないと動かない)
    double px = -sin(phi);
    double py = cos(phi);
    
    //TODO maxTime = nInfo.arraySize ??
    for(int i=0; i < maxTime; i++)
    {
      double complex WTH = 0 + 0 + Wz[k+i]*sz;
      double complex WPH = 0 + 0;
      double complex UTH = Ux[k+i]*sx + Uy[k+i]*sy + 0;
      double complex UPH = Ux[k+i]*px + Uy[k+i]*py;
      double complex ETH = coef*(-Z_0_S*WTH-UPH);
      double complex EPH = coef*(-Z_0_S*WPH+UTH);
      
      Eth[k+i] = ETH; //TODO : 物理単位に変換が必要かも
      Eph[k+i] = EPH;
    }
  }
}
      
//時間領域のEthの書き出し.
void ntffTM_TimeOutput(dcomplex *Ux, dcomplex *Uy, dcomplex *Wz)
{
  const int maxTime = field_getMaxTime();
  NTFFInfo nInfo = field_getNTFFInfo();
  dcomplex *Eth, *Eph;  
  Eth = newDComplex(360*nInfo.arraySize);
  Eph = newDComplex(360*nInfo.arraySize);

  //Eth, Ephを計算
  ntffTM_TimeTranslate(Ux,Uy,Wz,Eth,Eph);
  
  double **out_ref = (double**)malloc(sizeof(double*) * (LAMBDA_EN_NM-LAMBDA_ST_NM+1));

  for(int l=0; l<=LAMBDA_EN_NM-LAMBDA_ST_NM; l++)
    out_ref[l] = newDouble(360);

  //fft用に2の累乗の配列を確保
  dcomplex *eth = newDComplex(NTFF_NUM);
  for(int ang=0; ang<360; ang++)
  {
    int k= ang*nInfo.arraySize;

    memset((void*)eth,0,sizeof(dcomplex)*NTFF_NUM);               //0で初期化
    memcpy((void*)eth, (void*)&Eth[k], sizeof(dcomplex)*maxTime); //コピー
    cfft(eth, NTFF_NUM); //FFT

    FieldInfo fInfo = field_getFieldInfo();
    for(int lambda_nm=LAMBDA_ST_NM; lambda_nm<=LAMBDA_EN_NM; lambda_nm++)
    {
      //線形補完 TODO : index = n-1 となるほどの小さいlambdaを取得しようとするとエラー
      double p = C_0_S * fInfo.h_u_nm * NTFF_NUM / lambda_nm;
      int index = floor(p);
      p = p-index;
      out_ref[lambda_nm-LAMBDA_ST_NM][ang] = ((1-p)*cnorm(eth[index]) + p*cnorm(eth[index+1]))/NTFF_NUM;
    }
  }  
  freeDComplex(eth);
  
  //カレントディレクトリを取得
  char buf[256], parent[512];
  getcwd(parent, 512);

  // fft変換後のデータをテキストファイルで書き出し
  sprintf(buf, "%d[deg].txt",(int)field_getWaveAngle());
  ntff_outputEnormTxt(out_ref, buf);
  printf("saved %s/%s\n", parent, buf);

  // fft変換後のデータをバイナリファイルで書き出し
  sprintf(buf, "%d[deg]_%dnm_%dnm_b.dat", (int)field_getWaveAngle(), LAMBDA_ST_NM, LAMBDA_EN_NM);
  ntff_outputEnormBin(out_ref, buf);
  printf("saved %s/%s\n", parent, buf);

  //メモリの解放
  for(int l=0; l<=LAMBDA_EN_NM-LAMBDA_ST_NM;l++)
      freeDouble(out_ref[l]);
  free(out_ref);
  /*
  // fft変換前の時間サンプリングデータを書き出し
  char re[1024], im[1024];
  sprintf(re, "%d[deg]_Eth_r.txt", (int)field_getWaveAngle());
  sprintf(im, "%d[deg]_Eth_i.txt", (int)field_getWaveAngle());
  FILE *fpRe = openFile(re);
  FILE *fpIm = openFile(im);
  for(int ang=0; ang<360; ang++){
    int k= ang*nInfo.arraySize;
    for(int i=0; i < maxTime; i++){
      fprintf(fpRe,"%.20lf " , creal(Eth[k+i]));
      fprintf(fpIm,"%.20lf " , cimag(Eth[k+i]));
    }
    fprintf(fpRe,"\n");
    fprintf(fpIm,"\n");
  }
  printf("saved %s/ %s & %s \n", parent, re, im);
  fclose(fpRe);
  fclose(fpIm);
  */
  free(Eth);
  free(Eph);
}

// eとU[stp] もしくは hとW[stp]を渡す
//UW_ang = Ux[stp],Uy[stp],Wz[stp] の事. array[360][num]を一次元配列で表しており, その角度における配列を引数にとる
static inline void calc(double time_plus_timeShift, dcomplex eh,  dcomplex *UW_ang){  
  int m = floor(time_plus_timeShift+0.5);
  double a = (0.5 + time_plus_timeShift - m);
  double b = 1.0-a;
  double ab = a-b;
//  double coef = 1.0/(4*M_PI*C_0_S*R0);
  UW_ang[m-1] += eh*b;//*coef;
  UW_ang[m]   += eh*ab;//*coef;
  UW_ang[m+1] -= eh*a;//*coef;
}

//時間領域に置けるNTFF(領域分割なし)
//ここでは coef = 1.0/(4*M_PI*C_0_S*R0) を乗算していない.
//最後にかけた方が効率が良いので, 最後にかける必要がある.
void ntffTM_TimeCalc(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz)
{
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();  
//  const double coef = 1.0/(4*M_PI*C_0_S*R0);
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2
  
  NTFFInfo nInfo = field_getNTFFInfo();

  // 中心から各頂点へのベクトルを計算しておく
  double lt_cx = nInfo.left - nInfo.cx;
  double rt_cx = nInfo.right - nInfo.cx;
  double bm_cy = nInfo.bottom - nInfo.cy;
  double tp_cy = nInfo.top - nInfo.cy;

  //各頂点のインデックス(Hx,Hy,Ez用)も計算しておく
  int lb_index   = field_index(nInfo.left , nInfo.bottom);
  int lt_index   = field_index(nInfo.left , nInfo.top);
  int rb_index   = field_index(nInfo.right, nInfo.bottom);
  int rt_index   = field_index(nInfo.right, nInfo.top);
  
  double ToRad = M_PI/180.0;
  int index_ang = 0;//角度angの0番目のインデックス
  for(int ang=0; ang<360; ang++, index_ang+=nInfo.arraySize )
  {
    double rad = ang*ToRad;
    double r1x_per_c = cos(rad)/C_0_S, r1y_per_c = sin(rad)/C_0_S;

    //ang°の位置にシフトしたポジション, こうすれば Ux_ang[i]でその角度のi番目にアクセスできる.
    dcomplex *Ux_ang = &Ux[index_ang];
    dcomplex *Uy_ang = &Uy[index_ang];
    dcomplex *Wz_ang = &Wz[index_ang];
    
     // (l,b) -> (r,b)  W = Js = n × H = ( 0, 0, Hx)  U = Ms = E × n = (Ez, 0,  0)
    {
      const double r2x = lt_cx, r2y = bm_cy; // (lt-cx, bm-cy)
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC; // (dot(r1, r2) + Rf)/C
      for(int k=lb_index; k<rb_index; k+=fInfo_s.N_PY) //右端は次のループで計算するため < になる
      {
        calc(timeE+timeShift, Ez[k], Ux_ang);
        calc(timeH+timeShift, 0.5*(Hx[k]+Hx[k-1]), Wz_ang);
        timeShift -= r1x_per_c; //右に1セル進むと -r1x_per_c大きくなる(timeShiftの式より)
      }
    }
    
    // (r,b)->(r,t)  W = Js = n × H = (0, 0,Hy)  U = Ms = E × n = ( 0,Ez,0)
    {
      const double r2x = rt_cx;    double r2y = bm_cy;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC;
      for( int k=rb_index; k<rt_index; k++) {
        calc(timeE+timeShift, Ez[k], Uy_ang);
        calc(timeH+timeShift, 0.5 * ( Hy[k] + Hy[k-fInfo_s.N_PY] ), Wz_ang);
        timeShift -= r1y_per_c; //上に1セルあがると -r1y_per_c大きくなる
      }
    }

    //(l,t)->(r,t) normal(0,1,0)  //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    {
      double r2x  =  lt_cx;  const double r2y  = tp_cy;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC; // (dot(r1, r2) + Rf)/C
      for ( int k=lt_index; k<rt_index; k+=fInfo_s.N_PY) {
        calc(timeE+timeShift,               -Ez[k], Ux_ang);
        calc(timeH+timeShift, -0.5*(Hx[k]+Hx[k-1]), Wz_ang);
        timeShift -= r1x_per_c; //右に1セル進むと -r1x_per_c大きくなる
      }
    }

    // (l,b) -> (l,t)  normal(-1,0,0)   //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    {
      const double r2x = lt_cx;    double r2y = bm_cy;
      double timeShift = -(r1x_per_c*r2x + r1y_per_c*r2y) + nInfo.RFperC; // (dot(r1, r2) + Rf)/C
      for (int k=lb_index;  k<lt_index;  k++) {      
        calc(timeE+timeShift, -Ez[k], Uy_ang);
        calc(timeH+timeShift, -0.5*(Hy[k]+Hy[k-fInfo_s.N_PY] ), Wz_ang);
        timeShift -= r1y_per_c; //上に1セルあがると -r1y_per_c大きくなる
      }
    }
  }
}

/*
//==================================================
//以下分割領域バージョン

typedef struct SubNTFFInfo_S
{
  int SUB_L, SUB_R, SUB_T, SUB_B; //各辺の位置
  int SUB_CX, SUB_CY;
  //以下は各辺への1次元配列インデックス,範囲外は-1
  int rb_index, rt_index; //右辺
  int tl_index, tr_index; //上底
  int lb_index, lt_index; //左辺
  int bl_index, br_index; //下底
}SubNTFFInfo_S;

static IndexOfIntPath indexPath;

static NTFFInfo ntffInfo;

static SubNTFFInfo_S subNInfo_s;

static int between(int st, int en, int a)
{
  return st <= a && a <= en ? a : -1;
}

//コンストラクタ
SubNTFFInfo_S newSubNTFFInfo_S()
{
  SubNTFFInfo_S info;
  info.SUB_L = info.SUB_R = info.SUB_T = info.SUB_B = -1;
  info.rb_index = infor.rt_index = infor.tl_index = infor.tr_index = -1;
  info.lb_index = infor.lt_index = infor.bl_index = infor.br_index = -1;
  return info;
}

void ntff_init()
{
  FieldInfo_S fInfo_s = field_getFieldInfo_S();
  //領域から5セル内側の長方形を積分路とする
  ntffInfo.bottom = fInfo_s.N_PML + 5;
  ntffInfo.top    = fInfo_s.N_PY - fInfo_s.N_PML - 5;
  ntffInfo.left   = fInfo_s.N_PML + 5;
  ntffInfo.right  = fInfo_s.N_PX - fInfo_s.N_PML - 5;

  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    ntffInfo.top - subInfo_s.OFFSET_Y; //上面
  int bm = ntffInfo.bottom - subInfo_s.OFFSET_Y; //下面
  int rt = ntffInfo.right  - subInfo_s.OFFSET_X; //右
  int lt = ntffInfo.left   - subInfo_s.OFFSET_X; //左

  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  subNInfo_s = newSubNTFFInfo_S(); //初期化

  FieldInfo   fInfo   = field_getFieldInfo();
  subNInfo_s.SUB_CX = fInfo.N_PX/2 - subInfo_s.OFFSET_X;
  subNInfo_s.SUB_CY = fInfo.N_PY/2 - subInfo_s.OFFSET_X;
  // 左辺よりも左側の領域か, 右辺よりも右側の領域ならば, 積分路の完全外なのでそれの確認
  if(lt <= subInfo_s.SUB_N_PX-2 && rt >= 1)
  {
    subNInfo_s.SUB_T = between(1, subInfo_s.SUB_N_PY-2, tp);    
    subNInfo_s.SUB_B = between(1, subInfo_s.SUB_N_PY-2, bm);
    //1次元インデックスの計算
    if( subNInfo_s.SUB_B >= 0 ){
      subNInfo_s.bl_index = field_subIndex( max(1,lt), bm );
      subNInfo_s.br_index = field_subIndex( min(rt, subInfo_s.SUB_N_PX-2), bm);
    }
    if( subNInfo_s.SUB_T >=0 )    {
      subNInfo_s.tl_index = field_subIndex( max(1,lt), tp );
      subNInfo_s.tr_index = field_subIndex( min(rt, subInfo_s.SUB_N_PX-2), tp);
    }
  }
  
  if( bm <= subInfo_s.SUB_N_PY-2 && tp >= 1)
  {
    subNInfo_s.SUB_R = between(1, subInfo_s.SUB_N_PX-2, rt);
    subNInfo_s.SUB_L = between(1, subInfo_s.SUB_N_PX-2, lt);
    //1次元インデックスの計算
    if( subNInfo_s.SUB_N_L >= 0){
      subNInfo_s.lb_index = field_subIndex( lt, max(1,bm) );
      subNInfo_s.lt_index = field_subIndex( lt, min(subInfo_s.SUB_N_PY-2,tp) );
    }    
    if( subNInfo_s.SUB_N_R >= 0){
      subNInfo_s.lb_index = field_subIndex( rt, max(1,bm) );
      subNInfo_s.lt_index = field_subIndex( rt, min(subInfo_s.SUB_N_PY-2,tp) );
    }  
  }
}

void ntff_TMTime_MPI(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz)
{
  FieldInfo_S fInfo_s      = field_getFieldInfo_S();
  SubFieldInfo_S subInfo_s = field_getSubFieldInfo_S();
  WaveInfo_S waveInfo_s    = field_getWaveInfo_S();
  NTFFInfo nInfo = field_getNTFFInfo();
  
  const double coef = 1.0/(4*M_PI*C_0_S*R0);

  //分割領域系に置ける, 中心の位置
  double cx = fInfo_s.N_PX/2 - subInfo_s.OFFSET_X;
  double cy = fInfo_s.N_PY/2 - subInfo_s.OFFSET_Y;

  //分割領域のインデックスに変換 (0以下, SUB_N_PX, PY以上だと範囲外)
  int tp =    nInfo.top - subInfo_s.OFFSET_Y; //上面
  int bm = nInfo.bottom - subInfo_s.OFFSET_Y; //下面
  int rt = nInfo.right  - subInfo_s.OFFSET_X; //右
  int lt = nInfo.left   - subInfo_s.OFFSET_X; //左

  int m_e, m_h;
  double a_e, b_e, ab_e, a_h, b_h, ab_h;
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2

  for(int ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);
    int stp = ang*nInfo.arraySize;  //角度ごとのインデックス
    
    //bottom normal(0,-1,0)    //W = Js = n × H = ( 0, 0, Hx)  U = Ms = E × n = (Ez, 0,  0)
    if ( 0 < bm && bm < subInfo_s.SUB_N_PY-1 ) {
      //積分路の端で無ければ, 分割領域の端から端までが積分路である
      // 一番外側はのりしろなので 1~SUB_N_PX-2まで, ただし lt <= x < rt なので 1~SUB_N_PY-1になってる
      const int subLeft  = max(1, lt);
      const int subRight = min(subInfo_s.SUB_N_PX-1, rt);
      for(int i=subLeft; i<subRight; i++)
      {
        //原点との距離
        const double r2x = i  - cx;
        const double r2y = bm - cy;
        
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = field_subIndex(i, bm);
        const double complex ez = Ez[k];
        const double complex hx = 0.5 * ( Hx[k] + Hx[k-1] );
        
        Ux[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hx*b_h*coef;
        Ux[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hx*ab_h*coef;
        Ux[stp+m_e+1] -= ez*a_e*coef;
        Wz[stp+m_h+1] -= hx*a_h*coef;
      }
    }
    
    //right normal(1,0,0)    //Js = n × H = (0, 0,Hy)  Ms = E × n = ( 0,Ez,0)
    if ( 0 < rt && rt < subInfo_s.SUB_N_PX-1) {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ ) {
        const double r2x = rt-cx;
        const double r2y =  j-cy;

        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&rt, &j);
        const double complex ez = Ez[k];
        const double complex hy = 0.5 * ( Hy[k] + Hy[k-subInfo_s.SUB_N_PY] );
        
        Uy[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hy*b_h*coef;      
        Uy[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hy*ab_h*coef;      
        Uy[stp+m_e+1] -= ez*a_e*coef;              
        Wz[stp+m_h+1] -= hy*a_h*coef;
      }
    }

    //top normal(0,1,0)  //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    if ( 0 < tp && tp < subInfo_s.SUB_N_PY-1) {      
      int subLeft  = max(1, lt);
      int subRight = min(SUB_N_PX, rt);
      for ( int i=subLeft; i<subRight; i++ ) {
        const double r2x  =  i-cx;
        const double r2y  = tp-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC; 
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&i, &tp);
        const double complex ez = -Ez[k];
        const double complex hx = -0.5 * ( Hx[k] + Hx[k-1] );

        Ux[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hx*b_h*coef;
        Ux[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hx*ab_h*coef;
        Ux[stp+m_e+1] -= ez*a_e*coef;
        Wz[stp+m_h+1] -= hx*a_h*coef;
      }
    }

    // (left,top) -> (left,bottom)  normal(-1,0,0)   //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    if ( 0 < lt && lt < SUB_N_PX ) {
      int subTop    = min(SUB_N_PY, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ ) {
        const double r2x = lt-cx;
        const double r2y =  j-cy;
        const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
        ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
        ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

        int k = subInd(&lt, &j);
        const double complex ez = -Ez[k];
        const double complex hy = -0.5 * ( Hy[k] + Hy[k - subInfo_s.SUB_N_PY] );
      
        Uy[stp+m_e-1] += ez*b_e*coef;
        Wz[stp+m_h-1] += hy*b_h*coef;      
        Uy[stp+m_e]   += ez*ab_e*coef;
        Wz[stp+m_h]   += hy*ab_h*coef;      
        Uy[stp+m_e+1] -= ez*a_e*coef;      
        Wz[stp+m_h+1] -= hy*a_h*coef;
      }
    }
  }
}
*/
//周波数領域のNTFF(分割領域バージョン)
void ntffTM_FrequencySplit( dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex resultEz[360])
{
  FieldInfo fInfo      = field_getFieldInfo();
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

  for(int ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double rx  = cos(rad), ry = sin(rad); //方向ベクトル

    dcomplex Nz = 0;
    dcomplex Lx = 0;
    dcomplex Ly = 0;
		
    // (left,bottom) -> (right,bottom)
    // 法線ベクトルはn=(0, -1)
    //積分路が分割領域内にあるか確認 
    if (  0 < bm && bm < subInfo_s.SUB_N_PY-1)
    {
      int subLeft  = max(1, lt);                  //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subRight = min(subInfo_s.SUB_N_PX-1, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで

      for ( int i=subLeft; i<subRight; i++ )
      {
        double r2x  =  i-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        double r2y  = bm-cy;
        int k = field_subIndex(i, bm);
        dcomplex C_EZ = Ez[k];        
        dcomplex C_HX = 0.5 * ( Hx[k] + Hx[k-1] );
        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HX * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0   ,0,C_HX)
        Lx  += C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (C_EZ,0,   0)
      }
    }

    // (right,bottom) -> (right,top) n=(1,0)
    if ( 0 < rt && rt < subInfo_s.SUB_N_PX-1)
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm);
      for ( int j=subBottom; j<subTop; j++ )
      {
        double r2x  = rt-cx; //cx, cyを分割領域空間に変換してるので, そのまま引き算が出来る
        double r2y  =  j-cy;

        int k = field_subIndex(rt, j);
        dcomplex C_EZ = Ez[k];
        dcomplex C_HY = 0.5 * ( Hy[k] + Hy[k-subInfo_s.SUB_N_PY] );

        double innerProd = rx*r2x + ry*r2y;  //内積
        Nz  += C_HY * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,   0, C_HY)
        Ly  += C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (0,C_EZ,    0)
      }
    }


    // (right,top) -> (left,top)  n=(0,1)
    if ( 0 < tp && tp < subInfo_s.SUB_N_PY-1)
    {      
      //左右の端で無ければ, 分割領域の端から端までが積分路である
      int subLeft  = max(1, lt);
      int subRight = min(subInfo_s.SUB_N_PX-1, rt); //領域の端っこでなれば, 全部が積分路なのでSUB_N_PX-1 でなく SUB_N_PX まで      
      for ( int i=subLeft; i<subRight; i++ )
      {
        double r2x  =  i-cx;
        double r2y  = tp-cy;

        int k = field_subIndex(i, tp);
        dcomplex C_EZ = Ez[k];
        dcomplex C_HX = 0.5 * ( Hx[k] + Hx[k-1] );

        double innerProd = rx*r2x  + ry*r2y;  //内積
        Nz   -= C_HX * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,    0, -C_HX)
        Lx   -= C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (-C_EZ,0,     0)
      }
    }

    // (left,top) -> (left,bottom)
    if ( 0 < lt && lt < subInfo_s.SUB_N_PX-1 )
    {
      int subTop    = min(subInfo_s.SUB_N_PY-1, tp);
      int subBottom = max(1, bm); 
      for ( int j=subBottom; j<subTop; j++ )
      {
        double r2x = lt-cx;
        double r2y =  j-cy;

        int k = field_subIndex(lt, j);
        dcomplex C_EZ  = Ez[k];
        dcomplex C_HY  = 0.5 * ( Hy[k] + Hy[k - subInfo_s.SUB_N_PY] );

        double innerProd = rx*r2x  + ry*r2y;  //内積      
        Nz   -= C_HY * cexp( I * waveInfo_s.K_s * innerProd );  // J = n × H = (0,     0, -C_HY)		
        Ly   -= C_EZ * cexp( I * waveInfo_s.K_s * innerProd );  // M = E × n = (0, -C_EZ,     0)
      }
    }
    dcomplex Lphi = -Lx*sin(rad) + Ly*cos(rad); //極座標変換
    resultEz[ang] = coef * ( Z_0_S*Nz + Lphi ) * fInfo.h_u_nm;
  }
  
  FILE *fp = fopen("ntffStr.txt", "w");
  
  for(int ang = 0; ang<360; ang++)
  {
    fprintf(fp, "%.18lf\n", cnorm(resultEz[ang]));
  }
  printf("ntff tm frequency finished\n");
}


/*
  
static inline void ntffCoef(double time, double timeShift, int *m, double *a, double *b, double *ab)
{
  double t = time + timeShift;
  *m = floor(t + 0.5); //四捨五入
  *a = (0.5 + t - *m);
  *b = 1.0-*a;
  *ab = *a-*b;
}


void ntffTM_TimeCalc2(dcomplex *Hx, dcomplex *Hy, dcomplex *Ez, dcomplex *Ux, dcomplex *Uy, dcomplex *Wz)
{  
  NTFFInfo nInfo = field_getNTFFInfo();
  const double C = C_0_S;
  const double R = 1.0e6;
  const double coef = 1.0/(4*M_PI*C*R);

  //分割領域系に置ける, 中心の位置
  const int cx = N_PX/2;
  const int cy = N_PY/2;
  
  int m_e, m_h;
  double a_e, b_e, ab_e, a_h, b_h, ab_h;  
  double timeE = field_getTime() - 1;   //t - Δt
  double timeH = field_getTime() - 0.5; //t - Δt/2

  //分割領域における積分路, 有効範囲は(1~SUB_N_PX-2)まで
  int tp =    nInfo.top;
  int bm = nInfo.bottom;
  int rt = nInfo.right;
  int lt = nInfo.left;

  int ang;
  //360°方向の, 遠方界を求める
  for(ang=0; ang<360; ang++)
  {
    double rad = ang*M_PI/180.0;
    double r1x = cos(rad), r1y = sin(rad);
    const int stp = ang*nInfo.arraySize;  //角度ごとのインデックス
    dcomplex *Ux_ang = &Ux[stp];
    dcomplex *Uy_ang = &Uy[stp];
    dcomplex *Wz_ang = &Wz[stp];

    // (l,b) -> (r,b)  W = Js = n × H = ( 0, 0, Hx)  U = Ms = E × n = (Ez, 0,  0)
    for(int i=lt; i<rt; i++)
    {
      //原点との距離
      const double r2x = i  - cx;
      const double r2y = bm - cy;
        
      const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      int k = field_index(i, bm);
      const double complex ez = Ez[k];
      const double complex hx = 0.5 * ( Hx[k] + Hx[k-1] );
      calc(timeE+timeShift, ez, Ux_ang);
      calc(timeH+timeShift, hx, Wz_ang);
    }

    // (r,b)->(r,t)  W = Js = n × H = (0, 0,Hy)  U = Ms = E × n = ( 0,Ez,0)
    for ( int j=bm; j<tp; j++ )
    {
      const double r2x = rt-cx;
      const double r2y =  j-cy;

      const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      int k = field_index(rt, j);
      const double complex ez = Ez[k];
      const double complex hy = 0.5 * ( Hy[k] + Hy[k-N_PY] );
//      calc(timeE+timeShift, ez, Uy_ang);
//      calc(timeH+timeShift, hy, Wz_ang);
      
      Uy[stp+m_e-1] += ez*b_e*coef;
      Uy[stp+m_e]   += ez*ab_e*coef;
      Uy[stp+m_e+1] -= ez*a_e*coef;              
      Wz[stp+m_h-1] += hy*b_h*coef;      
      Wz[stp+m_h]   += hy*ab_h*coef;
      Wz[stp+m_h+1] -= hy*a_h*coef;      
    }

    //(l,t)->(r,t) normal(0,1,0)  //Js = n × H = (0, 0,-Hx)  Ms = E × n = (-Ez, 0,  0)
    for ( int i=lt; i<rt; i++ )
    {
      const double r2x  =  i-cx;
      const double r2y  = tp-cy;
      const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC; 
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      int k = field_index(i, tp);
      const double complex ez = -Ez[k];
      const double complex hx = -0.5 * ( Hx[k] + Hx[k-1] );
//      calc(timeE+timeShift, ez, Ux_ang);
//      calc(timeH+timeShift, hx, Wz_ang);
      
      Ux[stp+m_e-1] += ez*b_e*coef;
      Ux[stp+m_e]   += ez*ab_e*coef;
      Ux[stp+m_e+1] -= ez*a_e*coef;
      Wz[stp+m_h-1] += hx*b_h*coef;
      Wz[stp+m_h]   += hx*ab_h*coef;
      Wz[stp+m_h+1] -= hx*a_h*coef;
      
    }

    // (l,b) -> (l,t)  normal(-1,0,0)   //Js = n × H = (0,0,-Hy)    Ms = E × n = (0,-Ez,0)
    for ( int j=bm; j<tp; j++ )
    {
      const double r2x = lt-cx;
      const double r2y =  j-cy;
      const double timeShift = -(r1x*r2x + r1y*r2y)/C + nInfo.RFperC;
      ntffCoef(timeE, timeShift, &m_e, &a_e, &b_e, &ab_e);
      ntffCoef(timeH, timeShift, &m_h, &a_h, &b_h, &ab_h);

      int k = field_index(lt, j);
      const double complex ez = -Ez[k];
      const double complex hy = -0.5 * ( Hy[k] + Hy[k-N_PY] );
//      calc(timeE+timeShift, ez, Uy_ang);
//      calc(timeH+timeShift, hy, Wz_ang);
      
      Uy[stp+m_e-1] += ez*b_e*coef;
      Uy[stp+m_e]   += ez*ab_e*coef;
      Uy[stp+m_e+1] -= ez*a_e*coef;     
      Wz[stp+m_h-1] += hy*b_h*coef;     
      Wz[stp+m_h]   += hy*ab_h*coef;
      Wz[stp+m_h+1] -= hy*a_h*coef;
      
    }
  }
}

 */
