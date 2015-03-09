#ifndef STENCIL_TESTBENCH_H
#define STENCIL_TESTBENCH_H
#include<mpi.h>
#include "./stencil.h"

#define setpixel(x_,y_)           \
  {               \
    int xx,yy;              \
    xx=(x_)%sizex;            \
    yy=(y_)%sizey;            \
    xx=(xx+sizex)%sizex;          \
    yy=(yy+sizey)%sizey;          \
    v[yy*sizex+xx]=1;           \
  }

int main(int argc, char* argv[]);
void draw_circle(unsigned char *v, int sizex, int sizey, int x0, int y0, int r);
void write_pgm(FILE *f, unsigned char *v, int sizex, int sizey);

#endif
