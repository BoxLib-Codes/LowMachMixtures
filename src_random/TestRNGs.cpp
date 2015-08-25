#include <iostream>
#include <mpi.h>

extern "C" void SeedParallelRNG (int globalSeed);

int
main (int argc, char** argv)
{
    int NProcs, MyProc, IOProc = 0;
    
    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &NProcs);

    MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);

    if(MyProc == IOProc)
        std::cout << "Testing RNG seeding ..." << std::endl;

    SeedParallelRNG (1776);

    MPI_Finalize();
}
