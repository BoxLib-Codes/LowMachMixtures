#include <math.h>

#include "HydroGrid.h"
#include "RNGs.h"

#define NDIMS 3
#define NMAX 16

main ()
{
   int i,j,k;
   int Ncells = NMAX*NMAX*NMAX;
   double densities[NMAX*NMAX*NMAX];
   double velocities[NDIMS*NMAX*NMAX*NMAX];
   int nCells[3] = {NMAX,NMAX,NMAX}; // Two-dimensional grid
   double systemLength[3] = {NMAX*1.0, NMAX*1.0, NMAX*1.0}; // Cell size is 1.0
   double heatCapacity[1] = {1.0};
   double dt=1.0; // Time step between snapshots
   double std=0.01; // Standard deviation
   double old;
   
   srandgen (10426); // Seed the generator
   
   createHydroAnalysis_C (nCells, 1, NDIMS, 1, systemLength, heatCapacity, dt, 0, pow(std*std,-1), 1);

   for(j=0; j<Ncells; j++) {
      genrandn((densities+j));
      densities[j] = 1.0 + std*densities[j];
   }   
   
   for(i=0; i<1000; i++) // Number of timesteps
   {
      for(j=0; j<NDIMS*Ncells; j++) {
         genrandn((velocities+j));
         velocities[j]=std*velocities[j];
      }
      for(j=0; j<Ncells; j++) {
         old = densities[j];
         genrandn((densities+j));
         // We purposely put some memory here
         densities[j] = 1.0 + (0.5*std*densities[j] + 0.5*(old-1.0));
      }   
      updateHydroAnalysisIsothermal_C (velocities, densities);
   }
   
   writeToFiles_C(-1); // Write to files
   
   destroyHydroAnalysis_C ();

}
