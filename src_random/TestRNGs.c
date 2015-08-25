#include <mpi.h>

void SeedParallelRNG (int globalSeed); // In SeedRNGs.cpp

int
main (int argc, char** argv)
{
    int NProcs, MyProc, IOProc = 0;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);

    MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);

    SeedParallelRNG (1776);

    MPI_Finalize();
}
