#include "drawer.h" //ここに_USE_OPENGLを定義
#include <stdio.h>
#include <stdlib.h>
#ifdef _USE_OPENGL
  #include <GL/glew.h>
  #include <GLUT/glut.h>
#endif
#include <mpi.h>
#include <complex.h>
#include "simulator.h"
#include "mpiTM_UPML.h"
#include "mpiTE_UPML.h"
int windowX = 100;
int windowY = 100;
int windowWidth = 300;
int windowHeight=300;

void display(void)
{
#ifdef _USE_OPENGL
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
#endif
}

void idle(void)
{
  simulator_calc();

  if( simulator_isFinish() )
  {
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
    MPI_Finalize();
    exit(0);
  }
  
#ifdef _USE_OPENGL
  glutPostRedisplay();  //再描画
#endif  
  MPI_Barrier(MPI_COMM_WORLD);
}

int main( int argc, char *argv[] )
{
    MPI_Init( 0, 0 );
    
    int    width  = 2560; //横幅(nm)
    int    height = 2560; //縦幅(nm)
    double   h_u  = 10;   //1セルの大きさ(nm)
    int       pml = 10;  //pmlレイヤの数
    double lambda = 500;  //波長(nm)
    int      step = 2000; //計算ステップ
    enum MODEL   modelType = MIE_CYLINDER; // モデルの種類
    enum SOLVER solverType = TE_UPML_2D;        // 計算方法
    simulator_init(width, height, h_u, pml, lambda, step, modelType, solverType);    //simulator

#ifdef _USE_OPENGL    
    enum COLOR_MODE colorMode = CABS;
    
    glutInit(&argc, argv);
    glutInitWindowPosition(windowX,windowY);
    glutInitWindowSize(windowWidth, windowHeight);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE);
    glutCreateWindow("FDTD Simulator");
    glutDisplayFunc(display);
    glutIdleFunc(idle);
    glewInit();
    drawer_init(colorMode);
    glutMainLoop();
    MPI_Finalize();
#endif

#ifndef _USE_OPENGL    //only calculate mode

    while(!simulator_isFinish())
    {
       simulator_calc();    
    }
    MPI_Barrier(MPI_COMM_WORLD);
    simulator_finish();
    MPI_Finalize();
#endif

    return 1;
}
