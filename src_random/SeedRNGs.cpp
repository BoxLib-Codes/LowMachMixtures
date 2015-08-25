#include <iostream>
#include <cstdlib>
#include <set>
#include <vector>

#include <mpi.h>

extern "C"
{
#include "RNGs.h"
}

extern "C" void SeedParallelRNG (int globalSeed)
{
  int SEED;
  int NProcs, MyProc, IOProc = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &NProcs);

  MPI_Comm_rank(MPI_COMM_WORLD, &MyProc);

  std::vector<int> seeds(NProcs);

  SEED = globalSeed;
  
  if(MyProc == IOProc)
     std::cout << "Initial SEED = " << SEED << "\n\n";

  if (NProcs>1)
  {
    // This is based on Mike Lijewski's code in LLNS/main.cpp
    if(MyProc == IOProc)
    {
        //
        // We'll get NProcs unique integers with which to seed MT.
        //
        srand(SEED);

        std::set<int> seedset;

        do
        {
            SEED=rand();

            if (SEED)
                //
                // Don't consider zero a valid seed.
                //
                seedset.insert(SEED);
        }
        while (seedset.size() < NProcs);

        int i = 0;
        for (std::set<int>::const_iterator it = seedset.begin();
             it != seedset.end();
             ++it, ++i)
        {
            seeds[i] = *it;
        }

        for (int i = 0; i < NProcs; i++)
        {
            std::cout << "CPU# " << i << ", SEED = " << seeds[i] << '\n';
        }
        std::cout << '\n';

    }

    MPI_Scatter(&seeds[0],
                1,
                MPI_INT,
                &SEED,
                1,
                MPI_INT,
                IOProc,
                MPI_COMM_WORLD);
  }                
  //
  // Now set'm ...
  //
  srandgen(SEED);

}
