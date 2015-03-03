#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

void write_pgm(FILE *f, unsigned char *v, int sizex, int sizey)
{
  int i, j;

  fprintf(f, "P2\n");
  fprintf(f, "%d %d\n", sizex, sizey);
  fprintf(f, "%d\n", 255);

  for (i = 0; i < sizey; i++ )
  {
    for (j = 0; j < sizex; j++ )
    {
      fprintf(f, " %d", (unsigned int)v[i * sizex + j]);
    }
    fprintf(f, "\n");
  }
}


#define setpixel(x_,y_)           \
  {               \
    int xx,yy;              \
    xx=(x_)%sizex;            \
    yy=(y_)%sizey;            \
    xx=(xx+sizex)%sizex;          \
    yy=(yy+sizey)%sizey;          \
    v[yy*sizex+xx]=1;           \
  }


void smooth(unsigned char *v, int sizex, int sizey)
{
  int x, y;

  for ( x = 1; x < sizex - 1; x++ )
  {
    for ( y = 1; y < sizey - 1; y++ )
    {
      v[y * sizex + x] =
        ( 0.40 * v[y * sizex + x] +
          0.15 * v[(y - 1) * sizex + x] +
          0.15 * v[(y + 1) * sizex + x] +
          0.15 * v[y * sizex + (x + 1)] +
          0.15 * v[y * sizex + (x - 1)] );
    }
  }

}


void draw_circle(unsigned char *v, int sizex, int sizey,
                 int x0, int y0, int r)
{
  int f = 1 - r;
  int ddF_x = 1;
  int ddF_y = -2 * r;
  int x = 0;
  int y = r;


  //left, right, upper and lower bound for radius
  setpixel(x0 - r, y0);
  setpixel(x0 + r, y0);
  setpixel(x0  , y0 - r);
  setpixel(x0  , y0 + r);

  while (x < y)
  {
    if (f >= 0)
    {
      y--;
      ddF_y += 2;
      f += ddF_y;
    }
    x++;
    ddF_x += 2;
    f += ddF_x;
    setpixel(x0 + x, y0 + y);
    setpixel(x0 - x, y0 + y);
    setpixel(x0 + x, y0 - y);
    setpixel(x0 - x, y0 - y);
    setpixel(x0 + y, y0 + x);
    setpixel(x0 - y, y0 + x);
    setpixel(x0 + y, y0 - x);
    setpixel(x0 - y, y0 - x);
  }
}

int main(int argc, char* argv[])
{
  unsigned char* global_array = NULL;
  int sizex, sizey;
  int niter;
  FILE *outfile;

  sizex = 1000;
  sizey = 1000;
  niter = 20;

  MPI_Init(&argc, &argv);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (size != 16)
  {
    fprintf(stderr, "%s: Only works with np=%d for now\n", argv[0], 16);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  if (rank == 0)
  {
    global_array = malloc(sizex * sizey * sizeof(unsigned char));
    for (int i = 0; i < sizex * sizey; i++)
      global_array[i] = 255;

    outfile = fopen("testimg.pgm", "w");

    draw_circle(global_array, sizex, sizey, 0, 0, 40);
    draw_circle(global_array, sizex, sizey, 0, 0, 30);
    draw_circle(global_array, sizex, sizey, 100, 100, 10);
    draw_circle(global_array, sizex, sizey, 100, 100, 20);
    draw_circle(global_array, sizex, sizey, 100, 100, 30);
    draw_circle(global_array, sizex, sizey, 100, 100, 40);
    draw_circle(global_array, sizex, sizey, 100, 100, 50);
  }

  //split into 4x4 domain
  int ndims, lsizex, lsizey, nprocdim;
  nprocdim = 4;
  ndims = 2;
  lsizex = sizex / nprocdim;
  lsizey = sizey / nprocdim;

  int dims[2] = {nprocdim, nprocdim};
  int periods[2] = {0, 0};

  unsigned char* local_array = malloc(lsizex * lsizey * sizeof(unsigned char));

  for (int i = 0; i < nprocdim * nprocdim; i++)
  {
    local_array[i] = 0;
  }



  MPI_Comm cart_comm;
  int coords[2];
  int cart_rank;

  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, 0, &cart_comm);
  //MPI_Cart_coords(cart_comm, rank, ndims, coords);

  MPI_Datatype blocktype_temp, blocktype;
  int sizes[2]    = {sizex, sizey}; /* size of global array */
  int subsizes[2] = {lsizex, lsizey}; /* size of sub-region */
  int starts[2]   = {0, 0};  /* let's say we're looking at region "0",
                                 which begins at index [0,0] */

  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_UNSIGNED_CHAR, &blocktype_temp);
  MPI_Type_create_resized(blocktype_temp, 0, lsizex * sizeof(unsigned char), &blocktype);
  MPI_Type_commit(&blocktype);

  int sendcounts[nprocdim * nprocdim];
  int displs[nprocdim * nprocdim];

  if (rank == 0)
  {
    for (int i = 0; i < nprocdim * nprocdim; i++) sendcounts[i] = 1;
    int disp = 0;
    for (int i = 0; i < nprocdim; i++)
    {
      for (int j = 0; j < nprocdim; j++)
      {
        displs[i * nprocdim + j] = disp;
        disp += 1;
      }
      disp += (lsizex - 1) * nprocdim;
    }
  }

  // Send partitions out to the different processes
  MPI_Scatterv(global_array, sendcounts, displs, blocktype, local_array, lsizex * lsizey, MPI_UNSIGNED_CHAR, 0, cart_comm);

  /*
    //1. exchange guard cells
    int topLineStart = coords[1] * sizex + coords[0] * lsizex;
    int bottomLineStart = (coords[1] * lsizey + lsizey - 1);
    //int xend = xstart + block_size_x -1;

  */

  /*
  for (int i = 0; i < niter; i++ )
  {
    smooth(global_array, sizex, sizey);
  }

  */

  /* it all goes back to process 0 */
  MPI_Gatherv(local_array, lsizex * lsizey,  MPI_UNSIGNED_CHAR,
              global_array, sendcounts, displs, blocktype,
              0, cart_comm);



  free(local_array);

  if (rank == 0)
  {
    write_pgm(outfile, global_array, sizex, sizey);
    fclose(outfile);
    free(global_array);
  }

  MPI_Type_free(&blocktype);

  MPI_Finalize();

  return 0;
}
