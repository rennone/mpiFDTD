#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>
#include "simulator.h"
#include "models.h"
#include "field.h"
#include "parser.h"

typedef struct Config
{
  FieldInfo field_info;
  int startAngle, endAngle, deltaAngle;
  enum MODEL ModelType;
  enum SOLVER SolverType;
}Config;

#define ST_PHI -90
#define EN_PHI 0
#define DELTA_PHI 5

// 以下 OPEN_GLの関数
#ifdef USE_OPENGL

#include "drawer.h"
#define WINDOW_WIDTH 300
#define WINDOW_HEIGHT 300
#include <GL/glew.h>

//Macの場合
#ifdef MAC_OS
#include <GLUT/glut.h>
#endif

//Mac以外
#ifndef MAC_OS
#include <GL/glut.h>
#endif

//プロトタイプ宣言
static void drawField();
static void drawSubField();
static void display();
static void idle();

#endif

int rank;
int numProc;
Config config;
char root[512]; //ルートディレクトリ

//ファイルからのパラメータ呼び出し
static void readConfig(FieldInfo *field_info)
{
  FILE *fp = NULL;
  if( !(fp = fopen("config.txt", "r")) )
  {
    printf("cannot find config.txt of main\n");
  }

  int err;
  char buf[1024],tmp[1024];

  parser_nextLine(fp, buf);
  field_info->width_nm    = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->height_nm   = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->h_u_nm      = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->pml         = atoi(buf);

  parser_nextLine(fp, buf);
  field_info->lambda_nm   = strtod(buf,tmp);

  parser_nextLine(fp, buf);
  field_info->stepNum     = atoi(buf);

//入射角度
  parser_nextLine(fp, buf);
  config.startAngle = atoi(buf);

  parser_nextLine(fp, buf);
  config.endAngle = atoi(buf);

  parser_nextLine(fp, buf);
  config.deltaAngle = atoi(buf);

  //モデルの種類
  parser_nextLine(fp, buf);
  config.ModelType = atoi(buf);

  //Solver情報
  parser_nextLine(fp, buf);
  config.SolverType = atoi(buf);

  fclose(fp);
}

static void initConfigFromText()
{
  //同時にファイルを読み込むのはヤバい気がするので, rank0だけが読み込んで同期する
  if(rank == 0)
  {    
    readConfig(&config.field_info);
    printf("===========FieldSetting=======\n");
    printf("fieldSize(nm) = (%d, %d) \nh_u = %d \npml = %d\n", config.field_info.width_nm, config.field_info.height_nm,
           config.field_info.h_u_nm, config.field_info.pml);
    printf("lambda(nm) = %d  \nstep = %d\n", config.field_info.lambda_nm, config.field_info.stepNum);

    printf("angle = %d .. %d (delta = %d)\n", config.startAngle, config.endAngle, config.deltaAngle);
    printf("==============================\n");
  }
    
  //configを同期
  if(rank == 0)
  {
    for(int i=1; i<numProc; i++)
    {
      MPI_Send((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, i, 1, MPI_COMM_WORLD);
    }
  }else{
    MPI_Status status;
    MPI_Recv((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
  }
}

static void calcFieldSize(FieldInfo *fInfo, int x_nm, int y_nm)
{
  //モデルのサイズ + (pml + ntff)*2 + 余白
  fInfo->width_nm  = x_nm + fInfo->h_u_nm*(fInfo->pml + 5)*2 + 200;
  fInfo->height_nm = y_nm + fInfo->h_u_nm*(fInfo->pml + 5)*2 + 200;

  bool square = false;
  if(square)
  {
    if(fInfo->width_nm > fInfo->height_nm)
      fInfo->height_nm = fInfo->width_nm;
    else
      fInfo->width_nm = fInfo->height_nm;
  }
}

static void initParameter()
{
  config.field_info.h_u_nm    = 10;
  config.field_info.pml       = 10;
  config.field_info.lambda_nm = 500;
  config.field_info.stepNum   = 1500;
  config.startAngle = ST_PHI;
  config.endAngle   = EN_PHI;
  config.deltaAngle = DELTA_PHI;  
  config.SolverType = TM_UPML_2D;

  int x_nm, y_nm;
  models_needSize(&x_nm, &y_nm);
  calcFieldSize(&config.field_info, x_nm, y_nm);

  printf("%d, %d\n", config.field_info.width_nm, config.field_info.height_nm);
}

int main( int argc, char *argv[] )
{
  getcwd(root, 512); //カレントディレクトリを保存
  
  models_setModel(LAYER); //MIE_CYLINDER
  
  MPI_Init( 0, 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

//  initConfigFromText();
  initParameter();  //パラメータを設定
  
  //プロセスごとに角度を分ける
  config.field_info.angle_deg = config.startAngle + config.deltaAngle*rank;

  //必要以上の入射角度をしようとしてもスルー
  if(config.field_info.angle_deg > config.endAngle)
    {
      printf("rank%d is finish\n",rank);
      MPI_Finalize(); //プロセスごとにFinalizeしてもok
      exit(0); //call finalize before exit()
    }
  
  //シミュレーションの初期化.
  simulator_init(config.field_info, config.ModelType, config.SolverType);  

  //exitしたプロセスがあると停止してしまうのでMPI_Finalizeは使えない
//  MPI_Barrier(MPI_COMM_WORLD); //(情報表示がずれないように)全員一緒に始める
  printf("rank=%d, angle=%d\n",rank, config.field_info.angle_deg);
  
#ifndef USE_OPENGL
  //only calculate mode
  while(1)
  {
    while(1)
    {
      //シミュレーションをまわす
      while(!simulator_isFinish()) {
        simulator_calc();    
      }

      //シミュレーション終わったら, 入射角度を変えて再計算
      int angle = field_getWaveAngle()+config.deltaAngle*numProc;
      if( angle <= config.endAngle )
      {      
        simulator_reset();      
        field_setWaveAngle(angle);
      } else {
        simulator_finish();
        break;
      }
    }

    break;
    /*
    //一旦シミュレーションは終了する.
    simulator_finish();
    //モデルを変更して再計算する
    if ( !models_isFinish() ) {      
      moveDirectory(root);    //カレントディレクトリを元に戻す.
      models_moveDirectory(); //もう一度潜る
      int x_nm, y
      simulator_init(con); //モデルが変わったのでソルバーも再計算する.
    } else
    {
      break;
      }*/
  }
  MPI_Finalize(); //プロセスごとにFinalizeしてもok
#endif

  
#ifdef USE_OPENGL
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  
  int windowX = 1.0*subInfo.OFFSET_X / subInfo.SUB_N_PX * WINDOW_WIDTH;
  int windowY = 800-1.0*subInfo.OFFSET_Y/subInfo.SUB_N_PY * WINDOW_HEIGHT - WINDOW_HEIGHT;
  enum COLOR_MODE colorMode = CREAL;
  glutInit(&argc, argv);
  glutInitWindowPosition(windowX,windowY);
  glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
  glutCreateWindow("FDTD Simulator");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  glewInit();
  drawer_init(colorMode);
  glutMainLoop();
#endif

  return 1;
}



// OpenGL関連の実装部分
#ifdef USE_OPENGL

static void drawField()
{
  FieldInfo_S sInfo = field_getFieldInfo_S();
  drawer_paintImage(sInfo.N_PML, sInfo.N_PML, sInfo.N_PX-sInfo.N_PML, sInfo.N_PY-sInfo.N_PML, sInfo.N_PX, sInfo.N_PY,
                    simulator_getDrawingData());
  
  drawer_paintModel(sInfo.N_PML, sInfo.N_PML, sInfo.N_PX-sInfo.N_PML, sInfo.N_PY-sInfo.N_PML, sInfo.N_PX, sInfo.N_PY,
                    simulator_getEps());
}

static void drawSubField()
{
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  drawer_paintImage(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getDrawingData());
  drawer_paintModel(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getEps());  
}

static void display()
{
  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  drawField(); //分割しない全体領域の描画
  //drawSubField();  //分割された領域の描画
  
  drawer_draw();     
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
}

static void idle(void)
{
  simulator_calc();

  //シミュレーションが続くなら再描画して終わり
  if( !simulator_isFinish() ){
    glutPostRedisplay();  //再描画
    return;
  }

  //シミュレーションが終わった => 角度を変える
  int angle = field_getWaveAngle()+config.deltaAngle*numProc;
  if( angle <= config.endAngle )
  {      
    simulator_reset();
    field_setWaveAngle(angle);
    return;
  }

  //角度も終わったらシミュレーションは終わる
  simulator_finish();  
  
  //角度も終わった => モデルを変える
  if( models_isFinish() )
  {
    //モデルの変更もしない場合は終わる
    MPI_Finalize();
    exit(0);        
  }
  else
  {
    //モデルを変更して再計算する
    moveDirectory(root);    //カレントディレクトリを元に戻す.
    models_moveDirectory(); //もう一度潜る

    initParameter();  //パラメータを設定
  
    //プロセスごとに角度を分ける
    config.field_info.angle_deg = config.startAngle + config.deltaAngle*rank;

    //必要以上の入射角度をしようとしてもスルー
    if(config.field_info.angle_deg > config.endAngle)
      {
	printf("rank%d is finish\n",rank);
	MPI_Finalize(); //プロセスごとにFinalizeしてもok
	exit(0); //call finalize before exit()
      }  
    //シミュレーションの初期化.
    simulator_init(config.field_info, config.ModelType, config.SolverType);  
    }

}
// 以上 OPENGLの関数
#endif
