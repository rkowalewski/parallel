#include <stdlib.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define N 1000

void prefixsum_inplace(int* x, int n)
{

  int ithread;
  int nthreads;
  int* suma;
  #pragma omp parallel shared (suma) private(ithread,nthreads)
  {
    int i;
    int ithread = omp_get_thread_num();
    int nthreads = omp_get_num_threads();
    #pragma omp single
    {
      suma = (int*) malloc(n * sizeof(int));
      suma[0] = 0;
    }
    int sum = 0;
    #pragma omp for schedule(static)
    for (i = 0; i < n; i++)
    {
      sum += x[i];
      x[i] = sum;
    }
    suma[ithread + 1] = sum;
    #pragma omp barrier
    int offset = 0;
    for (i = 0; i < (ithread + 1); i++)
    {
      offset += suma[i];
    }
    #pragma omp for schedule(static)
    for (i = 0; i < n; i++)
    {
      x[i] += offset;
    }
  }
}

int main(int argc, char** argv)
{
  int x[N];
  int i;
  for (i = 0; i < N; i++)
  {
    x[i] = i;
  }

  prefixsum_inplace(x, N);

  for (i = 0; i < N; i++)
  {
    printf("%d ", x[i]);
  }
}
