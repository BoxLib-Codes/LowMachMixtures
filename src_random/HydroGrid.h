/*
A. Donev
C prototypes for interface to HydroGridModule.f90 via HydroGridCInterface.f90 
*/ 

 void setHydroInputFile_C(const char * filename);
 
 /* This will read input from the file hydroGridOptions.nml */
 void createHydroAnalysis_C (int nCells[3], int nSpecies, int nVelocityDimensions, \
     int isSingleFluid, double systemLength[3], double heatCapacity[], double timestep, \
     int nPassiveScalars, double structFactMultiplier, int project2D);
 
 void destroyHydroAnalysis_C ();

 void resetHydroAnalysis_C ();
 
 void updateHydroAnalysisConserved_C (double density[], double current[], double energy[]);

 void updateHydroAnalysisIsothermal_C (double velocity[], double density[]);

 void updateHydroAnalysisMixture_C (double velocity[], double density[], double concentration[]);

 /*mcai------start----------------*/
 void projectHydroGrid_C (double density[], double concentration[], char * filename, int id, int save_snapshot);
 void writeHydroGridMixture_C (double density[], double concentration[], char * filename, int id);
 /*mcai------end------------------*/

 /* For simplicity, we alias density as an advected scalar here, if grid%nPassiveScalars>0, otherwise ignore it
  * Note: You can pass NULL or a zero-sized array for vz or density if there is no data (it simply won't be touched) */
 void updateHydroAnalysisStaggered_C (int nGhost, int nPassiveScalars, double vx[], double vy[], double vz[], int nGhostScalars, double scalar[]);

 /* The id is an additional integer to append to file names, if positive */
 void writeToFiles_C(int id);

 /* This is actually in Tridiagonal.f90 */
 void SolveTridiagonal(double *a, double *b, double *c, double *r, double *u, int n, int periodic);
