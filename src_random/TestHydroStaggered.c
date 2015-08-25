#include <stdlib.h>
#include <math.h>

#include "HydroGrid.h"
#include "RNGs.h"

#define NDIMS 2 // Two or three dimensional code?
#define NX 16 // Grid size (x)
#define NY 32 // Grid size (y)
#define NSCALS 1 // Zero if there is only velocities

main ()
{
   int i,j,k;
   int Ncells = NX*NY;
   double scalars[NX*NY*NSCALS];
   double vx[(NX+1)*NY];
   double vy[NX*(NY+1)];
   int nCells[3] = {NX,NY,1}; // Two-dimensional grid
   double systemLength[3] = {NX*1.0, NY*1.0, 1.0}; // Cell size is 1.0
   double heatCapacity[1] = {1.0};
   double dt=1.0; // Time step between snapshots
   double std=0.01; // Standard deviation of noise
   
   srandgen (10426); // Seed the generator
   
   createHydroAnalysis_C (nCells, 1, NDIMS, 1, systemLength, heatCapacity, dt, 0, pow(std*std,-1), 0);
   
   for(i=0; i<1000; i++) // Number of timesteps
   {
      for(j=0; j<(NX+1)*NY; j++) {
         genrandn((vx+j));
         vx[j]=std*vx[j];
      }
      for(j=0; j<NX*(NY+1); j++) {
         genrandn((vy+j));
         vy[j]=std*vy[j];
      }

      for(j=0; j<Ncells*NSCALS; j++) {
         genrandn((scalars+j));
         scalars[j] = 1.0 + std*scalars[j];
      }   
      
      updateHydroAnalysisStaggered_C (1, NSCALS, vx, vy, NULL, 0, scalars);
   }
   
   writeToFiles_C(-1); // Write to files
   
   destroyHydroAnalysis_C ();

}
