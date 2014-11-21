#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define MYTIMEVAL( tv_ )      \
  ((tv_.tv_sec)+(tv_.tv_usec)*1.0e-6)

#define MYTIMESTAMP( time_ ) \
{                             \
  static struct timeval tv;   \
  gettimeofday( &tv, NULL );  \
  time_=MYTIMEVAL(tv);        \
}

#define NUM_ITERATIONS 1000000000
#define NUM_DIGITS 10

int main(int argc, char *argv[])
{
  long i;
  double pi = 0;
  double tstart, tstop, time;
  //printf ("Max Threads: %d\n", omp_get_max_threads());
  //printf ("Num Threads: %d\n", omp_get_num_threads());

  MYTIMESTAMP(tstart);

  #pragma omp parallel for reduction(+: pi)
  for (i = 0; i < NUM_ITERATIONS; i++)
  {
    pi += 1.0 / (i * 4.0 + 1.0);
    pi -= 1.0 / (i * 4.0 + 3.0);
  }
  pi = pi * 4.0;

  MYTIMESTAMP(tstop);
  time = tstop - tstart;
  printf("Pi = %.*f\n", NUM_DIGITS, pi);
  printf("Time: %.*f\n", 2, tstop - tstart);
  return 0;
}
