#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>
#include "simulator.h"
#include "field.h"
#include "parser.h"

typedef struct Config
{
  int startAngle, endAngle, deltaAngle;
  enum MODEL ModelType;
  enum SOLVER SolverType;
}Config;

int rank;
int numProc;
Config config;

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

static void drawField()
{
  FieldInfo_S sInfo = field_getFieldInfo_S();
  drawer_paintImage(0,0, sInfo.N_X, sInfo.N_Y, sInfo.N_PX, sInfo.N_PY,
                    simulator_getDrawingData());
  drawer_paintModel(0,0, sInfo.N_X, sInfo.N_Y, sInfo.N_PX, sInfo.N_PY,
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
  
  if( simulator_isFinish())
  {
    MPI_Barrier(MPI_COMM_WORLD);
    if( field_getWaveAngle() < config.endAngle )
    {      
      simulator_reset();
      int angle = field_getWaveAngle()+config.deltaAngle*numProc;
      field_setWaveAngle(angle);
    }else{
      simulator_finish();
      MPI_Finalize();
      exit(0);
    }
  }
  glutPostRedisplay();  //再描画
//  MPI_Barrier(MPI_COMM_WORLD);
}
// 以上 OPENGLの関数
#endif

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

int main( int argc, char *argv[] )
{
  MPI_Init( 0, 0 );
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);

  FieldInfo field_info;
  //同時にファイルを読み込むのはヤバい気がするので, rank0だけが読み込んで同期する
  if(rank == 0)
  {    
    readConfig(&field_info);
    
    printf("===========FieldSetting=======\n");
    printf("fieldSize(nm) = (%d, %d) \nh_u = %d \npml = %d\n", field_info.width_nm, field_info.height_nm,
           field_info.h_u_nm, field_info.pml);
    printf("lambda(nm) = %d  \nstep = %d\n", field_info.lambda_nm, field_info.stepNum);

    printf("angle = %d .. %d (delta = %d)\n", config.startAngle, config.endAngle, config.deltaAngle);
    printf("==============================\n");
  }

  //configを同期
  if(rank == 0)
  {
    for(int i=1; i<numProc; i++)
    {
      MPI_Send((int*)&field_info, sizeof(FieldInfo)/sizeof(int), MPI_INT, i, 0, MPI_COMM_WORLD);
      MPI_Send((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, i, 1, MPI_COMM_WORLD);
    }
  }else{
    MPI_Status status;
    MPI_Recv((int*)&field_info, sizeof(FieldInfo)/sizeof(int), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    MPI_Recv((int*)&config, sizeof(Config)/sizeof(int), MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
  }
  
  field_info.angle_deg = config.startAngle + config.deltaAngle*rank;

  //シミュレーションの初期化. 同時に, 必要なディレクトリまで移動している.
  simulator_init(field_info, config.ModelType, config.SolverType);  
  
  MPI_Barrier(MPI_COMM_WORLD); //(情報表示がずれないように)全員一緒に始める
  printf("rank=%d, angle=%d\n",rank, field_info.angle_deg);
#ifdef USE_OPENGL
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  
  int windowX = 1.0*subInfo.OFFSET_X / subInfo.SUB_N_PX * WINDOW_WIDTH;
  int windowY = 800-1.0*subInfo.OFFSET_Y/subInfo.SUB_N_PY * WINDOW_HEIGHT - WINDOW_HEIGHT;
  enum COLOR_MODE colorMode = CABS;
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
//  MPI_Finalize();
#endif

#ifndef USE_OPENGL
  //only calculate mode
  while(1)
  {
    while(!simulator_isFinish()) {
      simulator_calc();    
    }
    
//    MPI_Barrier(MPI_COMM_WORLD);
    if( field_getWaveAngle() < config.endAngle )
    {      
      simulator_reset();
      int angle = field_getWaveAngle()+config.deltaAngle*numProc;
      field_setWaveAngle(angle);
    } else {
      simulator_finish();
      break;
    }
  }
  MPI_Finalize(); //プロセスごとにFinalizeしてもok
#endif

  return 1;
}
