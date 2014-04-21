#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <complex.h>
#include "simulator.h"

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

void display(void)
{

  glEnableClientState( GL_VERTEX_ARRAY );
  glEnableClientState( GL_TEXTURE_COORD_ARRAY );

  int Nx, Ny, Npx, Npy;
  simulator_getSubFieldPositions(&Nx, &Ny, &Npx,&Npy );
  drawer_paintImage(1,1, Nx, Ny, Npx, Npy, simulator_getDrawingData());
  drawer_paintModel(1,1, Nx, Ny, Npx, Npy, simulator_getEps());
  
  drawer_draw();
    
  glDisableClientState( GL_VERTEX_ARRAY );
  glDisableClientState( GL_TEXTURE_COORD_ARRAY );
  glutSwapBuffers();
}

void idle(void)
{
  simulator_calc();

  if( simulator_isFinish() )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
//    MPI_Finalize();
    exit(0);
  }
  glutPostRedisplay();  //再描画
//  MPI_Barrier(MPI_COMM_WORLD);
}
#endif

int main( int argc, char *argv[] )
{
  MPI_Init( 0, 0 );
  FieldInfo field_info;
  field_info.width_nm  = 2560;
  field_info.height_nm = 2560;
  field_info.h_u_nm    = 10;
  field_info.pml       = 10;
  field_info.lambda_nm = 500;
  field_info.stepNum   = 2000;
  enum MODEL   modelType = MIE_CYLINDER; // モデルの種類
  enum SOLVER solverType = MPI_TM_UPML_2D;        // 計算方法
  simulator_init(field_info, modelType, solverType);

#ifdef USE_OPENGL
  int windowX = 100;
  int windowY = 100;
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
