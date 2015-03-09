#ifndef STENCIL_H
#define STENCIL_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define TAG_EXCHANGE_GHOSTS 9

void calc_stencil(unsigned char* array, int sizex, int sizey, int niter, int ndims);
void smooth(unsigned char *src, unsigned char *dst, int x, int y, int sizex);

#endif
