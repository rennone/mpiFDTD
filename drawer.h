#ifndef _DRAWER_H
#define _DRAWER_H
#include "myComplex.h"

enum COLOR_MODE
{
  CREAL,
  CABS
};

#ifdef USE_OPENGL
extern void (*drawer_getDraw(void))(void);
extern void drawer_paintImage(int l, int b, int r, int t,int wid, int hei, double complex*);
extern void drawer_paintModel(int l, int b, int r, int t,int wid, int hei, double *);
extern void drawer_paintTest(void);
extern void drawer_init(enum COLOR_MODE);
extern void drawer_finish(void);
extern void drawer_draw(void);
#endif //USE_OPENGL

extern void drawer_outputImage(char *fileName, dcomplex* data, double *model, int width, int height);
extern void drawer_screenshot(const char* filename);

#endif //_DRAWER_H
