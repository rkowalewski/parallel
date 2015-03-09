#include "./stencil.h"

void calc_stencil(unsigned char *array, int sizex, int sizey, int nprocs, int rank, int ndims)
{
  //split into 4x4 domain
  int by, bx, px, py;
  //Get Dimensions
  int dims[2] = {0, 0};
  MPI_Dims_create(nprocs, ndims, dims);

  px = dims[0];
  py = dims[1];

  if (!rank)
  {
    if ((sizey % py != 0) || (sizex % px != 0))
    {
      fprintf(stderr, "size not evenly divided by number of processors per dimension\n");
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
  }

  by = sizey / py;
  bx = sizex / px;


  //Create Cartesian Topology
  MPI_Comm cart_comm;
  int reorder = 0;
  int periods[2] = {0, 0};
  MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &cart_comm);

  int coords[2];
  MPI_Cart_coords(cart_comm, rank, ndims, coords);

  int west, north, east, south;
  //north to south
  MPI_Cart_shift(cart_comm, 1, 1, &north, &south);
  //west to east
  MPI_Cart_shift(cart_comm, 0, 1, &west, &east);

  MPI_Datatype array_scatter_type;
  int sizes[2]    = {sizey, sizex}; /* size of global array */
  int subsizes[2] = {by, bx}; /* size of sub-region */
  int starts[2] = {0, 0};

  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_UNSIGNED_CHAR, &array_scatter_type);
  MPI_Type_create_resized(array_scatter_type, 0, sizeof(unsigned char), &array_scatter_type);
  MPI_Type_commit(&array_scatter_type);

  unsigned char *subarray, *aold, *anew, *tmp;
  subarray = malloc(by * bx * sizeof(unsigned char));
  aold = calloc((by + 2) * (bx + 2), sizeof(unsigned char));
  anew = calloc((by + 2) * (bx + 2), sizeof(unsigned char));


  int displs[nprocs];
  int sendcounts[nprocs];

  if (rank == 0)
  {
    for (int r = 0; r < nprocs; r++)
    {
      sendcounts[r] = 1;
      int proc_coords[2];
      MPI_Cart_coords(cart_comm, r, ndims, proc_coords);
      displs[r] = proc_coords[1] * px * bx * by + proc_coords[0] * bx;
    }
  }

  MPI_Scatterv(array, sendcounts, displs, array_scatter_type, subarray, bx * by, MPI_UNSIGNED_CHAR, 0, cart_comm);
  //put subarray in array with ghost cells
  for (int row = 1 ; row < by + 1; row++)
  {
    for (int col = 1; col < bx + 1; col++)
    {
      aold[row * (bx + 2) + col] = subarray[(row - 1) * bx + (col - 1)];
    }
  }

  // create east-west type
  MPI_Datatype north_south_type;
  MPI_Type_contiguous(bx, MPI_UNSIGNED_CHAR, &north_south_type);
  MPI_Type_commit(&north_south_type);
  // create east-west type
  MPI_Datatype west_east_type;
  MPI_Type_vector(by, 1, bx + 2, MPI_UNSIGNED_CHAR, &west_east_type);
  MPI_Type_commit(&west_east_type);

  MPI_Request reqs[8];

#define ind(i,j) (j) * (bx + 2) + (i)
  int i, j, iter;

  int northRow, southRow, westCol, eastCol;
  northRow = southRow = westCol = eastCol = -1;
  int *northRowPtr, *southRowPtr, *westColPtr, *eastColPtr;
  northRowPtr = &northRow;
  southRowPtr = &southRow;
  eastColPtr = &eastCol;
  westColPtr = &westCol;
  if (north != MPI_PROC_NULL)
  {
    *northRowPtr = 1;
  }
  else
  {
    northRowPtr = southRowPtr;
  }

  if (south != MPI_PROC_NULL)
  {
    *southRowPtr = by;
  }
  else
  {
    southRowPtr = northRowPtr;
  }

  if (west != MPI_PROC_NULL)
  {
    *westColPtr = 1;
  }
  else
  {
    westColPtr = eastColPtr;
  }

  if (east != MPI_PROC_NULL)
  {
    *eastColPtr = bx;
  }
  else
  {
    eastColPtr = westColPtr;
  }

  for (iter = 0; iter < NITER; iter++)
  {
    //start halo communication
    MPI_Isend(&aold[ind(1, by)] /* south */, 1, north_south_type, south, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[0]);
    MPI_Isend(&aold[ind(1, 1)] /* north */, 1, north_south_type, north, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[1]);
    MPI_Isend(&aold[ind(bx, 1)] /* east */, 1, west_east_type, east, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[2]);
    MPI_Isend(&aold[ind(1, 1)] /* west */, 1, west_east_type, west, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[3]);
    MPI_Irecv(&aold[ind(1, 0)] /* north */, 1, north_south_type, north, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[4]);
    MPI_Irecv(&aold[ind(1, by + 1)] /* south */, 1, north_south_type, south, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[5]);
    MPI_Irecv(&aold[ind(0, 1)] /* east */, 1, west_east_type, west, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[6]);
    MPI_Irecv(&aold[ind(bx + 1, 1)] /* west */, 1, west_east_type, east, TAG_EXCHANGE_GHOSTS, cart_comm, &reqs[7]);

    //update inner cells
    for ( int i = 2; i < by; i++)
    {
      for ( int j = 2; j < bx; j++ )
      {
        smooth(aold, anew, j, i, bx + 2);
      }
    }
    //wait for halo
    MPI_Waitall(8, reqs, MPI_STATUSES_IGNORE);

    //update outer cells (north - south)
    if (*northRowPtr > -1 && *southRowPtr > -1)
    {
      for (i = *northRowPtr; i <= *southRowPtr; i += by - 1)
      {
        for (j = 2; j < bx; j++)
        {
          smooth(aold, anew, j, i, bx + 2);
        }
      }
    }

    //update outer cells (east - west)
    if (*westColPtr > -1 && *eastColPtr > -1)
    {
      for (i = 1; i < by + 1; i++)
      {
        for (j = *westColPtr; j <= *eastColPtr; j += by - 1)
        {
          smooth(aold, anew, j, i, bx + 2);
        }
      }
    }

    if (iter < NITER - 1)
    {
      //swap arrays
      tmp = anew;
      anew = aold;
      aold = tmp;
    }
  }

  //Create array type to send back the smoothed subarrays without ghost cells
  MPI_Datatype array_gather_type;
  sizes[0] = by + 2;
  sizes[1] = bx + 2;
  subsizes[0] = by;
  subsizes[1] = bx;
  starts[0] = 1;
  starts[1] = 1;
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_UNSIGNED_CHAR, &array_gather_type);
  MPI_Type_create_resized(array_gather_type, 0, sizeof(unsigned char), &array_gather_type);
  MPI_Type_commit(&array_gather_type);

  MPI_Gatherv(anew, 1,  array_gather_type,
              array, sendcounts, displs, array_scatter_type,
              0, cart_comm);

  free(subarray);
  free(aold);
  free(anew);
  MPI_Type_free(&array_scatter_type);
  MPI_Type_free(&array_gather_type);
  MPI_Type_free(&north_south_type);
  MPI_Type_free(&west_east_type);
}

void smooth(unsigned char *src, unsigned char *dst, int x, int y, int sizex)
{
  dst[y * sizex + x] =
    ( 0.40 * src[y * sizex + x] +
      0.15 * src[(y - 1) * sizex + x] +
      0.15 * src[(y + 1) * sizex + x] +
      0.15 * src[y * sizex + (x + 1)] +
      0.15 * src[y * sizex + (x - 1)] );
}
