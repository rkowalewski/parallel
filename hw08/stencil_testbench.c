#include "./stencil_testbench.h"

int main(int argc, char* argv[])
{
  unsigned char* array = NULL;
  int sizex, sizey;
  int niter;
  FILE *outfile;

  MPI_Init(NULL, NULL);

  int rank, size;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  if (argc < 3)
  {
    if (!rank) printf("usage: stencil <n> <niters>\n");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  if (!rank)
  {
    sizex = sizey = atoi(argv[1]);
    niter = atoi(argv[2]);
    int args[2] = {sizex, niter};

    array = malloc(sizex * sizey * sizeof(unsigned char));

    for (int i = 0; i < sizex * sizey; i++)
      array[i] = 255;

    outfile = fopen("testimg.pgm", "w");

    draw_circle(array, sizex, sizey, 0, 0, 40);
    draw_circle(array, sizex, sizey, 0, 0, 30);
    draw_circle(array, sizex, sizey, 100, 100, 10);
    draw_circle(array, sizex, sizey, 100, 100, 20);
    draw_circle(array, sizex, sizey, 100, 100, 30);
    draw_circle(array, sizex, sizey, 100, 100, 40);
    draw_circle(array, sizex, sizey, 100, 100, 50);

    MPI_Bcast(args, 2, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
  {
    int args[2];
    MPI_Bcast(args, 2, MPI_INT, 0, MPI_COMM_WORLD);
    sizex = sizey = args[0];
    niter = args[1];
  }


  double t = -MPI_Wtime();
  calc_stencil(array, sizex, sizey, size, rank, 2);
  t += MPI_Wtime();

  if (!rank)
  {
    printf("processing this image took %1.2f seconds\n", t);
    write_pgm(outfile, array, sizex, sizey);
    fclose(outfile);
    free(array);
  }


  MPI_Finalize();

  return EXIT_SUCCESS;
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
