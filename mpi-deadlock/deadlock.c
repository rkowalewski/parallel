#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>


void illustrate_MPI_send(char op, char* message, char mode, int world_rank);

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    printf("usage: deadlock [SRA] msglen [SBRV]");
    return EXIT_FAILURE;
  }

  char* op = argv[1];
  int len = atoi(argv[2]);
  char* mode = argv[3];

  char* message = malloc(len + 1);

  for (int idx = 0; idx < len - 1; idx++)
  {
    message[idx] = 'x';
  }

  message[len - 1] = '\0';

  MPI_Init(NULL, NULL);

  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // We are assuming at least 2 processes for this task
  if (world_size != 2)
  {
    fprintf(stderr, "World size must be two for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  illustrate_MPI_send(op[0], message, mode[0], world_rank);

  MPI_Finalize();

  free(message);

  return EXIT_SUCCESS;
}

void illustrate_MPI_send(char op, char* message, char mode, int world_rank)
{
  if (op == 'S' && mode == 'S')
  {
    if (world_rank == 0) {
    
    }
  }
}
