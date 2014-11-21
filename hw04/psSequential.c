#include <stdio.h>
#define N 1000
int main(int argc, char* argv[])
{
  int i;
  int A[N];
  int B[N];
  for (i = 0; i < N; i++)
  {
    A[i] = i;
  }
  B[0] = A[0];

  for (i = 1; i < N; i++)
  {
    B[i] = B[i - 1] + A[i];
  }
  for (i = 0; i < N; i++)
  {
    printf("%d ", B[i]);
  }
}
