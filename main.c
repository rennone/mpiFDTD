#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>
#include "simulator.h"
#include "field.h"
#include "parser.h"

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

void drawField()
{
  FieldInfo_S sInfo = field_getFieldInfo_S();
  drawer_paintImage(0,0, sInfo.N_X, sInfo.N_Y, sInfo.N_PX, sInfo.N_PY,
                    simulator_getDrawingData());
  drawer_paintModel(0,0, sInfo.N_X, sInfo.N_Y, sInfo.N_PX, sInfo.N_PY,
                    simulator_getEps());
}

void drawSubField()
{
  SubFieldInfo_S subInfo = field_getSubFieldInfo_S();
  drawer_paintImage(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getDrawingData());
  drawer_paintModel(1,1, subInfo.SUB_N_X, subInfo.SUB_N_Y, subInfo.SUB_N_PX, subInfo.SUB_N_PY,
                    simulator_getEps());  
}

void display()
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

void idle(void)
{
  simulator_calc();

  if( simulator_isFinish() )  {
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
//    MPI_Finalize();
    exit(0);
  }
  glutPostRedisplay();  //再描画
//  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

//ファイルからのパラメータ呼び出し
void readField(FILE *fp, FieldInfo *field_info)
{  
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
}

int main( int argc, char *argv[] )
{
  FILE *fp = NULL;
  if( !(fp = fopen("config.txt", "r")) )
  {
    printf("cannot find config.txt of main\n");
  }

  FieldInfo field_info;
  readField( fp, &field_info);

  MPI_Init( 0, 0 );
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  char buf[1024], tmp[1024];
  enum MODEL  modelType = MORPHO_SCALE;//MIE_CYLINDER; // モデルの種類  
  enum SOLVER solverType;

  printf("===========FieldSetting=======\n");
  printf("fieldSize(nm) = (%d, %d) \nh_u = %d \npml = %d\n", field_info.width_nm, field_info.height_nm,
         field_info.h_u_nm, field_info.pml);
  printf("lambda(nm) = %lf  \nstep = %d\n", field_info.lambda_nm, field_info.stepNum);
  printf("==============================\n");
  if(rank == 0)
    solverType  = TM_UPML_2D;        // 計算方法
  else
    solverType  = TE_UPML_2D;        // 計算方法
  
  simulator_init(field_info, modelType, solverType);

  
  field_setWaveAngle(-90);
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
  while(!simulator_isFinish())
  {
    simulator_calc();    
  }
//  MPI_Barrier(MPI_COMM_WORLD);
  simulator_finish();
//  MPI_Finalize();
#endif

  return 1;
}
