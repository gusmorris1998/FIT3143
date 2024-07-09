#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <time.h>
#include "shared.h"

int periods[DIMENSIONS];
int dims[DIMENSIONS];
int n, m, baseRank;
int size, cartSize;

int main(int argc, char **argv)
{
    int rank, provided;
    int coord[DIMENSIONS];

    // Initialisingh the comms
    MPI_Comm tempComm, node_grid_comm;

    // Initialization of thread 
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int periods[DIMENSIONS] = {0, 0}; // Non-periodic grid

    baseRank = cartSize = size - 1;

    // Split the comms into designated colour
    MPI_Comm_split(MPI_COMM_WORLD, rank == baseRank, 0, &tempComm);

    // Read Arguments
    dims[0] = atoi(argv[1]);
    dims[1] = atoi(argv[2]);

    char *dirname;
    dirname = "outputFolder";

    if (rank == baseRank) {
        runBase(dirname);
    }
    else {
        initialiseNodes(rank, tempComm, dirname);
        runNodes(MPI_COMM_WORLD, dirname);
        MPI_Comm_free(&tempComm);
    }

    MPI_Finalize();
    return 0;
}