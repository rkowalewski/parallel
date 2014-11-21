#include <stdio.h>
#include <sys/time.h>
#include <time.h>

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

#define NUM_ITERATIONS 1000000

int main(int argc, char *argv[])
{
  long i;
  long count = 0;
  double tstart, tstop, time;
  //printf ("Max Threads: %d\n", omp_get_max_threads());
  //printf ("Num Threads: %d\n", omp_get_num_threads());

  MYTIMESTAMP(tstart);
  
  #pragma omp parallel for reduction(+: count)
  for (i = 2; i <= NUM_ITERATIONS; i++) 
  {
    count += isprime(i);
  }

  MYTIMESTAMP(tstop);
  time = tstop-tstart;
  printf("prime count = %ld\n", count);
  printf("Time: %.*f\n", 2, time);
  return 0;
}

int isprime(int p)
{
  int d;
  for (d = 2; d < p; d++)
  {
    if (p % d == 0)
      return 0;
  }
  return 1;
}
