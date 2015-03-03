
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
  unsigned char* array;
  int i, sizex, sizey;
  int niter;
  FILE *outfile;

  sizex = 1000;
  sizey = 1000;
  niter = 20;

  array = malloc(sizex * sizey * sizeof(unsigned char));

  for (i = 0; i < sizex * sizey; i++)
    array[i] = 255;

  outfile = fopen("testimg.pgm", "w");

  draw_circle(array, sizex, sizey, 0, 0, 40);
  draw_circle(array, sizex, sizey, 0, 0, 30);
  draw_circle(array, sizex, sizey, 100, 100, 10);
  draw_circle(array, sizex, sizey, 100, 100, 20);
  draw_circle(array, sizex, sizey, 100, 100, 30);
  draw_circle(array, sizex, sizey, 100, 100, 40);
  draw_circle(array, sizex, sizey, 100, 100, 50);

  
  for ( i = 0; i < niter; i++ )
  {
    smooth(array, sizex, sizey);
  } 

  write_pgm(outfile, array, sizex, sizey);

  fclose(outfile);

  return 0;
}
