#ifndef STENCIL_H
#define STENCIL_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#define NITER 20
#define TAG_EXCHANGE_GHOSTS 9

void calc_stencil(unsigned char* array, int sizex, int sizey, int nprocs, int rank, int ndims);
void smooth(unsigned char *src, unsigned char *dst, int x, int y, int sizex);

#endif
