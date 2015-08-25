!----------------------------------------------------------------------------------------------------------------------------------!

!-------------------------------------------------------------
! Code for calculating statistical steady-state correlations (covariance matrix)
! of hydrodynamic values defined on a rectangular grid of shape (Nx,Ny,Nz)
! Ergodicity is assumed, that is, this code calculates time averages instead of ensemble averages
! Both real-space cross-correlations and Fourier structure factors are calculated for selected pairs of variables
! Written by Anton de la Fuente (SJSU) and Aleksandar Donev (Courant Institute), starting 2010
!-------------------------------------------------------------

!-------------------------------------------------------------
! For now we assume an ideal mixture of ideal gases, but it should not be hard to add other EOSs
 
! It is essentially assumed here that the fluctuations are equilibrium, or at least that
! there are implicit delta(r1-r2) and delta(k1+k2) in the covariance matrix of the hydrodynamic fluctuations
! do that only correlations between variables in the same cell or the same wavenumber are calculated
! This restriction implies that non-homogeneous states such as systems with walls cannot be fully analyzed

!-------------------------------------------------------------
! To be consistent with BoxLib, the hydrodynamic arrays have a shape (Nx,Ny,Nz,1:nVariables,0:nComponents-1)
! The component, or fluid, indexed with zero here always denotes the whole mixture (e.g., total momentum density over all species)
! nComponents can either be nSpecies when each species has a separate velocity/temperature (multi-fluid model), or
! nComponents=1 if there is a common velocity/temperature (single-fluid model)
! To be consistent with LLNS codes, the hydrodynamic variables are ordered in the following way:
! density, velocity (current), temperature (energy), concentration, passive scalars

!-------------------------------------------------------------
! A variable, such as rho_1, is specified as a (species,variable) pair, and species zero is the total
! So rho is (species=0,variable=1), but vy_1 is (species=1,variable=3)

! Single-fluid models:
!----------------------------------------
! If isSingleFluid=T, then nComponents=1 and the hydrodynamic variables are stored in an array of dimensions
! (Nx,Ny,Nz, nVariables = 1+nDimensions+1+(nSpecies-1)+nPassiveScalars, 0)
! Primitive variables are the total density, current and energy, and concentrations
! Conserved variables are densities of all the species, the total current, and total energy
! We explicitly store the partial densities since we want unbiased means of concentration

! For example, a single-species two-dimensional fluid with one passive scalar s would have the variables:
!     Primitive=(rho,vx,vy,T,s), Conserved=(rho,jx,jy,e,s)
! A two-species three-dimensional fluid with no passive scalars would have:
!     Primitive=(rho,vx,vy,vz,T,c1=rho1/rho), Conserved=(rho,jx,jy,jz,e,rho1)

! Multi-fluid models (no passive scalars)
!----------------------------------------
! If isSingleFluid=F, then nComponents=nSpecies and the hydrodynamic variables are stored in an array of dimensions
! (Nx,Ny,Nz, nVariables = 1+nDimensions+1+1, 0:nSpecies-1)
! For a single-species fluid there is no difference with the single-fluid models.
! For multi-species mixtures, primitive variables are the density, velocity, temperature and concentration of each species
! The concentration for species=0 is defined to be concentration of the last species
! Conserved variables are densities, current and energy for each of the species

! For example, a two-fluid mixture in two dimensions would have the variables:
!     Primitive: species=0: (rho,vx,vy,T,c2), species=1: (rho1,vx1,vy1,T1,c1)
!     Conserved: species=0: (rho,jx,jy,e), species=1: (rho1,jx1,jy1,e1)

!-------------------------------------------------------------
! TODOS:
! 1. Make a destruction routine for the grid, which will call a destructor for the structure factors
! 2. In the output file for the structure factor one should down in the header which variables the column refers to
! 3. Save fields in higher dimension as VTK files, instead of just one-dimensional profiles
! 4. Add routine to write out dynamic structure factor and improve the format of saving static factors 
! 5. Also calculate and output instantaneous means using the spatially-averaged instantaneous values, to reduce 1/N bias

module HydroGridModule
   use FFTW
   use VisIt_Writer
   use ISO_C_BINDING
   implicit none
   save
   integer, parameter, public :: nMaxDims = 3 ! One can treat 2D or 1D as a degenerate case
   integer, parameter, private :: wp = kind (0.0d0) ! Working precision (typically double)
   integer, parameter, private :: r_fft = wp
   
   real(wp), parameter, private :: pi = 3.1415926535897932384626433832795028842_wp
   
   integer, parameter, private :: nMaxCharacters = 1024 ! Dimension of character strings in declarations   
   
   type HydroGrid
      integer :: nDimensions = nMaxDims ! How many velocities are stored, 1 (vx), 2 (vx+by) or 3 (vx+vy+vz)
         ! Note that the actual grid is always treated as 3D, since 1D and 2D are simply degenerate cases with Nz/Ny=1
         
      ! Fluid variables
      ! -------------------------------
      ! What variables are actually represented for the fluid can vary:
      integer :: nSpecies = 1, nPassiveScalars = 0
      logical :: isSingleFluid = .true.
         ! Does each of the species have separate velocities or temperatures
         ! or is the whole fluid considered together + concentrations         
      integer :: nFluids = 0 ! nFluids=nComponents-1
         ! For single-fluid models nFluids=0, otherwise nFluids=nSpecies-1

      real (wp), allocatable :: heatCapacity(:) ! Of dimension (1:nSpecies)
         ! By passive scalars here we simply mean a quantity that does not affect the thermodynamics
         ! If other constitutive laws are added in the future this should be kept consistent...
      
      ! Shortcuts for easy indexing of arrays:
      integer :: mIdx = 1 ! Mass density index
      integer :: jIdx1 = 2 ! Current start index
      integer :: jIdx2 = 4 ! Current end index
      integer :: eIdx = 5 ! Energy index
      integer :: cIdx1 = 6 ! Concentration start
      integer :: cIdx2 = 5 ! Concentration end 
      integer :: sIdx1 = 6 ! Passive scalars start
      integer :: sIdx2 = 5 ! Passive scalars end 
      integer :: nVariables=5 ! No concentration or passive scalars by default: Note nVariables=sIdx2 (redundant)
      integer :: nConservedVariables=5 ! For actual "conserved" variables (in the absence of collisions)

      ! System dimensions:
      ! -------------------------------
      integer :: nCells(nMaxDims)=1 ! Grid dimensions (single level rectangular grid only)
      real (wp) :: systemLength(nMaxDims)

      ! File I/O:
      ! -------------------------------
      character (nMaxCharacters) :: outputFolder = "./"
      character (nMaxCharacters) :: filePrefix = "HydroAnalysis"

      ! Boundary conditions      
      ! -------------------------------
      logical :: periodic=.true. ! Is the special direction (axisToPrint) periodic or not
      integer :: axisToPrint = 2
      real(wp) :: topConcentration=0.5_wp, bottomConcentration=0.5_wp

      ! Mean values:
      ! -------------------------------
      logical :: writeMeansVTK=.false., writeSnapshotVTK=.false. ! Write to VTK      
      logical :: storeConserved = .true. ! Keep track of the mean conserved variables as well     
      integer :: staggeredVelocities = 0 ! Are the velocities on a staggered grid by +/-dx/2 (+1 or -1)?
      integer :: iStep = 0 ! Number of elapsed time steps (snapshots)

      ! We always internally store the primitive values, however, we always
      ! keep track of the mean conserved variables
      ! These are all of dimension (Nx,Ny,Nz,nVariables,0:nFluids)
      real (wp), allocatable :: primitive(:, :, :, :, :) 
      real (wp), allocatable :: meanPrimitive(:, :, :, :, :)
      real (wp), allocatable :: meanConserved(:, :, :, :, :)
      
      ! Variance (real-space pair correlations)
      ! -------------------------------
      ! We could keep track of a number of real-space cross-correlations between specified pairs of variables
      logical :: writeVariancesVTK=.false. ! Write to VTK      
      integer :: nVariances = 0 ! How many pairs of variables to calculate real-space cross-correlations for
      integer :: correlationCell(nMaxDims) = 1 ! If calculating correlations between cells, this cell is the reference one
      integer, allocatable :: variancePairs(:, :, :) ! Of dimension (2,2,nVariances)
         ! The first dimension stores the variable as a (species,variable) pair
         ! If the first variable is negative then calculates correlations with the reference cell (correlationCell)
      real (wp), allocatable :: variances(:, :, :, :) ! Dimension (Nx,Ny,Nz,nVariances)
         ! The actual covariances for pairs of variables, in the same cell, or between a cell and the reference cell
      
      ! Static structure factors (spatial spectral pair correlations)
      ! -------------------------------

      ! Internal work-space for FFTW library:
      integer(fftw_p) :: staticFFTplan=0 ! FFTW handle
      complex (r_fft), allocatable :: staticFFTarray(:, :, :) ! Dimension (Nx,Ny,Nz)
      
      ! The user can select which values to calculate FFTs of:
      integer :: nVariablesToFFT = 0 ! How many FFTs actually need to be calculated
      integer, allocatable :: variablesToFFT(:, :) ! Dimension (2,nVariablesToFFT)
         ! The first dimension stores the variable as a (species,variable) pair
      
      ! We need to store these so we don't have to recalculate the same FFT if it appears in multiple pairs of variables:
      complex (r_fft), allocatable :: staticFourierTransforms(:, :, :, :) ! Dimension (Nx,Ny,Nz,nVariablesToFFT)
         ! This is an internal temporary array in order to avoid 

      ! Cross-correlations between structure factors for different variables: 
      logical :: writeSpectrumVTK=.false. ! Write to VTK
      logical :: estimateCovariances=.false. ! Estimate the covariance of all pairs of spectra
      logical :: subtractMeanFT=.true. ! Subtract the mean of the structure factors (for non-uniform systems)
      real (wp) :: structFactMultiplier = 1.0_wp ! Normalization for structure factors
      integer :: nStructureFactors = 0 ! How many pairs of variables to calculate Fourier-space cross-correlations for
      integer, allocatable :: vectorFactors(:) ! Dimension (nStructureFactors), 1=diagonal, -1=off-diagonal
         ! Whether this component is part of the vector (velocity) structure factor (which is itself a tensor)
         ! >0 for diagonal component, <0 for off-diagonal (numbered along diagonals), 0 for non-velocity
      complex (r_fft), allocatable :: meanStaticTransforms(:, :, :, :) ! Dimension (Nx,Ny,Nz,nVariablesToFFT)
         ! The mean of the spatial FT for each of the needed variables, <SpatialFT(A)>
      integer, allocatable :: structureFactorPairs(:, :) ! Dimension (2,nStructureFactors)
         ! Which pairs of variables to calculate cross-structure factors for
         ! A variable here is no longer a (species,variable) pair, instead,
         ! it is an index (1:nVariablesToFFT) into staticFourierTransforms
      complex (r_fft), allocatable :: structureFactors(:, :, :, :) ! Dimension (Nx,Ny,Nz,nStructureFactors)
         ! Time-averaged structure factors, i.e., cross-correlation between the spectrum of pairs of variables A and B:
         ! <SpatialFT(A)*CONJ(SpatialFT(B))> (at the *same* wavenumber k)
         ! What is actually output is the covariance <SpatialFT(A)*CONJ(SpatialFT(B)>-<SpatialFT(A)>*CONJ(<SpatialFT(B)>)
      real (r_fft), allocatable :: structureFactorsCovariance(:, :, :, :)
         ! Dimension (Nx,Ny,Nz,nStructureFactors*(nStructureFactors+1)/2)         

      ! Dynamic structure factors (spatio-temporal spectral pair correlations)
      ! For each pair of variables for which we calculated static structure factors,
      ! we can calculate dynamic ones at selected wavenumbers
      ! -------------------------------      
      integer :: nSavedSnapshots = 1000 ! How many temporal snapshots to save per timeseries for each wavenumber
      integer :: iSample = 0 ! Current temporal snapshot (sample) >=0 and <=nSavedSnapshots
      integer :: iTimeSeries = 0 ! How many timeseries we have had so far
      real (wp) :: timestep = 1.0_wp ! Time between snapshots

      ! FFTW internal data
      integer(fftw_p) :: dynamicFFTplan=0 ! FFTW handle
      complex (r_fft), allocatable :: dynamicFFTarray(:) ! Dimension (nSavedSnapshots)
      
      ! The shape of the discrete k (wavenumber) grid:
      integer :: minWavenumber(nMaxDims)
      integer :: maxWavenumber(nMaxDims)
      
      ! Dynamic structure factors for selected k's and all omega's
      integer :: nWavenumbers = 0 ! How many wavenumbers to keep track (record) of
      integer :: correlationWavenumber = -1 ! If positive, this calculates temporal correlations
         ! Between each of the tracked wavenumbers and the selected wavenumber
      integer :: minWavefrequency, maxWavefrequency
      integer, allocatable :: selectedWavenumbers(:, :) ! Dimension (2/3,nWavenumbers): Indexes in k-grid of selected wavenumbers
      integer, allocatable :: selectedWaveindices(:, :) ! Conversion to FFT indices
      complex (r_fft), allocatable :: savedStructureFactors(:, :, :) ! Dimension (nSavedSnapshots,nWavenumbers,nVariablesToFFT)
         ! Saved history of the static factors for the selected wavenumbers
      complex (r_fft), allocatable :: dynamicFactors(:, :, :)
         ! Dimension (nSavedSnapshots,nWavenumbers,nStructureFactors)
         ! Time-averaged dynamic structure factors, <SpatioTemporalFT(A)*CONJ(SpatioTemporalFT(B))>
         ! What is actually output is the covariance:
         ! <SpatioTemporalFT(A)*CONJ(SpatioTemporalFT(B))> - <SpatioTemporalFT(A)>*CONJ(<SpatioTemporalFT(B)>)
      complex (r_fft), allocatable :: meanDynamicFactors(:, :, :) ! Dimension (nSavedSnapshots,nWavenumbers,nVariablesToFFT)
         ! The mean of the spatio-temporal FT for each of the needed variables, <SpatioTemporalFT(A)>

   end type HydroGrid
   
   logical, public, save :: writeAbsValue=.true. ! Only write absolute value of static structure factors to VTK files
   integer, public, save :: writeTheory=-1 ! Write the theoretical prediction for incompressible hydro (-2=for MD analysis, -1=none, 0=continuum, 1=MAC)
   integer, public, save :: useRelativeVelocity=-1 ! Calculate v12 to test two-fluid models
   
contains

!------------------------------------------------------------------------
! Creation and Initialization
!------------------------------------------------------------------------
subroutine createHydroAnalysis (grid, nCells, nSpecies, nVelocityDimensions, isSingleFluid, &
      systemLength, heatCapacity, timestep, &
      nPassiveScalars, fileUnit, fileNamePrefix, structFactMultiplier)
   type (HydroGrid), target :: grid
   integer, intent(in) :: nCells(nMaxDims), nSpecies, nVelocityDimensions
   real (wp), intent(in) :: systemLength(nMaxDims), timestep
   real (wp), intent(in), optional :: heatCapacity(nSpecies)
   logical, intent(in) :: isSingleFluid
   integer, intent(in), optional :: nPassiveScalars
   integer, intent(in), optional :: fileUnit ! If positive, read the namelist from the given file unit
   character(len=*), intent(in), optional :: fileNamePrefix ! Over-rides the one from the namelist if present!
   real (wp), intent(in), optional :: structFactMultiplier

   integer :: nameListFile, staggeredVelocities
   logical :: storeConserved

   integer :: nVariances, correlationWavenumber
   integer, dimension(nMaxDims) :: correlationCell
   character (nMaxCharacters) :: variancePairs

   integer :: nStructureFactors
   character (nMaxCharacters) :: structureFactorPairs, vectorStructureFactor

   integer :: nWavenumbers
   character (32*nMaxCharacters) :: selectedWavenumbers ! There may be lots of these
   integer :: nSavedSnapshots

   character (nMaxCharacters) :: outputFolder
   character (nMaxCharacters) :: filePrefix
   integer :: axisToPrint
   logical :: periodic

   logical :: writeSpectrumVTK, writeVariancesVTK, writeMeansVTK, writeSnapshotVTK, &
      estimateCovariances, subtractMeanFT

   integer :: iDimension, iWavenumber
   
   real :: topConcentration, bottomConcentration

   ! User input controls:      
   !---------------------------
   nameListFile = -1      
   if(present(fileUnit)) nameListFile=fileUnit
   if(nameListFile<=0) then   
      nameListFile = 114
      open (nameListFile, file = "hydroGridOptions.nml", status="old", action="read")
      call readHydroGridNML(nameListFile)
      close (nameListFile)                     
   else
      call readHydroGridNML(nameListFile)
   end if            

   ! Basics: 
   !---------------------------      
   grid%nCells = nCells            
   grid%nSpecies = nSpecies
   grid%isSingleFluid = isSingleFluid
   grid%nDimensions = nVelocityDimensions
   grid%storeConserved = storeConserved
   grid%staggeredVelocities = staggeredVelocities
   if(present(nPassiveScalars)) grid%nPassiveScalars=nPassiveScalars
   write(6,*)" "
   write(*,*) "Initializing HydroAnalysis for grid of size ", grid%nCells, &
      " nSpecies/isSingleFluid=", grid%nSpecies, grid%isSingleFluid

   grid%systemLength = systemLength
   grid%timestep = timestep

   allocate( grid%heatCapacity(nSpecies) )
   if(present(heatCapacity)) then
      grid%heatCapacity = heatCapacity
   else
      grid%heatCapacity = 1.0_wp
   end if   
   
   grid%estimateCovariances = estimateCovariances
   grid%subtractMeanFT = subtractMeanFT

   ! BCs
   grid%periodic = periodic
   grid%axisToPrint = axisToPrint
   grid%topConcentration = topConcentration
   grid%bottomConcentration = bottomConcentration

   ! File I/O:
   grid%writeSpectrumVTK=writeSpectrumVTK
   grid%writeVariancesVTK=writeVariancesVTK
   grid%writeMeansVTK=writeMeansVTK
   grid%writeSnapshotVTK=writeSnapshotVTK
   grid%outputFolder = outputFolder
   grid%filePrefix = filePrefix      
   if(present(fileNamePrefix)) grid%filePrefix=fileNamePrefix ! Override the namelist input
   
   if(present(structFactMultiplier)) grid%structFactMultiplier=structFactMultiplier

   ! The variable indexes depend on dimensionality etc:
   !---------------------------      
   grid%mIdx  = 1
   grid%jIdx1 = grid%mIdx + 1
   grid%jIdx2 = grid%jIdx1 + grid%nDimensions - 1
   grid%eIdx  = grid%jIdx2 + 1      

   ! If there is only one species, cIdx2=eIdx and so there is no concentration variables
   ! If the mixture is treated as one fluid then we need nSpecies-1 concentrations
   ! Otherwise we also store concentrations as primitive variables, cIdx2=eIdx+1
   ! The concentration for species=0 is defined to be concentration of the last species
   grid%cIdx1 = grid%eIdx + 1
   if (grid%isSingleFluid) then
      grid%cIdx2 = grid%eIdx + grid%nSpecies - 1
      grid%nFluids = 0
   else if(grid%nSpecies>1) then ! Calculate concentrations
      grid%cIdx2 = grid%eIdx + 1
      grid%nFluids = grid%nSpecies - 1
   else ! Single species mixture, no need to store concentrations
      grid%cIdx2 = grid%eIdx
      grid%nFluids = 0
   end if

   ! For the single-fluid formulation we can add some passive scalars:
   grid%sIdx1 = grid%cIdx2 + 1
   if (grid%isSingleFluid) then
      grid%sIdx2 = grid%sIdx1 + grid%nPassiveScalars - 1
   else
      grid%sIdx2 = grid%sIdx1 - 1     
      grid%nPassiveScalars = 0 ! It does not make sense to do passive scalars for the multi-fluid formulation
   end if   

   if (grid%isSingleFluid) then
      grid%nVariables = grid%nSpecies + grid%nDimensions + 1 + grid%nPassiveScalars
   else
      grid%nVariables = 1 + grid%nDimensions + 1 + 1
   end if   

   write(*,*) grid%nVariables, " Primitive variables indexes as: rho=", grid%mIdx, &
      " v=", grid%jIdx1,"-",grid%jIdx2, " T=", grid%eIdx, &
      " c=", grid%cIdx1,"-",grid%cIdx2, " s=", grid%sIdx1,"-",grid%sIdx2
      
   if (grid%isSingleFluid) then
      ! Conserved variables are densities of all the species, the total current, and total energy
      ! We explicitly store the partial densities since we want unbiased means of concentration
      grid%nConservedVariables = grid%nSpecies + grid%nDimensions + 1
   else
      ! Conserved variables are densities, current and energy for each of the species
      grid%nConservedVariables = 1 + grid%nDimensions + 1
   end if            

   ! Allocate storage for means:
   !---------------------------      
   allocate( grid%primitive     (grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nVariables, 0:grid%nFluids) )
   allocate( grid%meanPrimitive (grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nVariables, 0:grid%nFluids) )
   if(grid%storeConserved) then
      allocate( grid%meanConserved (grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nConservedVariables, 0:grid%nFluids) )
   end if   

   ! Real-space correlations:
   !---------------------------      
   grid%nVariances = nVariances
   if(any(correlationCell<=0)) then ! Choose the central cell by default
      correlationCell = (nCells+1)/2
   end if
   grid%correlationCell = correlationCell
   allocate(grid%variances(grid%nCells(1), grid%nCells(2), grid%nCells(3), nVariances))
   allocate(grid%variancePairs(2, 2, nVariances)) ! First dimension gives (species,variable) pair
   read(variancePairs, *) grid%variancePairs

   ! Fourier-space correlations:
   !---------------------------      
   grid%nStructureFactors = nStructureFactors
   grid%nWavenumbers = nWavenumbers
   grid%nSavedSnapshots = nSavedSnapshots
   grid%correlationWavenumber = correlationWavenumber

   call createStaticFactors(grid, structureFactorPairs, vectorStructureFactor)
   call createDynamicFactors(grid, selectedWavenumbers)

   ! Initialize various counters:
   !---------------------------      
   call resetHydroAnalysis(grid)      

contains   

   subroutine readHydroGridNML(nameListFile)
      integer, intent(inout) :: nameListFile

      namelist / hydroAnalysisOptions / storeConserved, variancePairs, nVariances, correlationCell, &
         nStructureFactors, structureFactorPairs, vectorStructureFactor, writeAbsValue, &
         nWavenumbers, selectedWavenumbers, nSavedSnapshots, axisToPrint, staggeredVelocities, &
         outputFolder, filePrefix, writeSpectrumVTK, writeVariancesVTK, writeMeansVTK, writeSnapshotVTK, &
         estimateCovariances, topConcentration, bottomConcentration, useRelativeVelocity, &
         writeTheory, periodic, subtractMeanFT, correlationWavenumber

      ! Setup default values:
      storeConserved=.true.
      variancePairs=""
      nVariances=0
      correlationCell=0 ! This will default to the central cell
      correlationWavenumber=-1 ! Not used unless set positive
      nStructureFactors=0
      structureFactorPairs=""
      vectorStructureFactor=""
      nWavenumbers=0
      selectedWavenumbers=""
      nSavedSnapshots=grid%nSavedSnapshots
      axisToPrint=grid%axisToPrint
      periodic=grid%periodic
      outputFolder=grid%outputFolder
      filePrefix=grid%filePrefix
      writeSpectrumVTK=.false.
      writeVariancesVTK=.false.
      writeMeansVTK=.false.
      writeSnapshotVTK=.false.
      estimateCovariances=.false.
      subtractMeanFT=.true.
      staggeredVelocities=grid%staggeredVelocities
      topConcentration=0.5_wp
      bottomConcentration=0.5_wp

      ! Now read the user inputs:   
      read (nameListFile, nml = hydroAnalysisOptions)

   end subroutine   

end subroutine createHydroAnalysis

subroutine createStaticFactors(grid, structureFactorPairsString, vectorStructureFactor)
   type (HydroGrid), target :: grid
   character (*), intent(in) :: structureFactorPairsString, vectorStructureFactor

   integer :: variablesToFFT(2, 2*grid%nStructureFactors)
   integer :: iPair, iStructureFactor, iVariable, jVariable, nVariablesToFFT
   integer :: structureFactorPairs(2, 2, grid%nStructureFactors)
   integer :: nFactors, iDimension, iWavenumber, iTemp

   nFactors = grid%nStructureFactors
   
   if (nFactors <= 0) return

   allocate(grid%structureFactorPairs(2, nFactors))
   allocate(grid%vectorFactors(nFactors))
   allocate(grid%structureFactors(grid%nCells(1), grid%nCells(2), grid%nCells(3), nFactors))

   read(structureFactorPairsString, *) structureFactorPairs
   if(vectorStructureFactor /= "") then
      read(vectorStructureFactor, *) grid%vectorFactors
   else
      grid%vectorFactors = 0
   end if   

   iVariable = 1
   do iStructureFactor = 1, nFactors
   do iPair = 1, 2
      variablesToFFT(:, iVariable) = structureFactorPairs(:, iPair, iStructureFactor)
      iVariable = iVariable + 1
   end do
   end do

   nVariablesToFFT = 2*nFactors
   do iVariable = 1, 2*nFactors
      do jVariable = iVariable + 1, 2*nFactors
         if (any(variablesToFFT(:, jVariable) < 0)) cycle
         if (all(variablesToFFT(:, iVariable) == variablesToFFT(:, jVariable))) then
            variablesToFFT(:, jVariable) = -1
            nVariablesToFFT = nVariablesToFFT - 1
         end if
      end do
   end do
   grid%nVariablesToFFT = nVariablesToFFT

   allocate(grid%variablesToFFT(2, nVariablesToFFT))
   allocate(grid%staticFourierTransforms(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nVariablesToFFT))
   ! Because we often do some special manipulations to the structure factors, it is better to store their mean explicitly:
   if(.true.) allocate(grid%meanStaticTransforms(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nVariablesToFFT))
   if(grid%estimateCovariances) then      
      allocate(grid%structureFactorsCovariance(grid%nCells(1), grid%nCells(2), grid%nCells(3), nFactors*(nFactors+1)/2))
   end if   

   jVariable = 1
   do iVariable = 1, 2*nFactors
      if (any(variablesToFFT(:,iVariable) < 0)) cycle
      grid%variablesToFFT(:, jVariable) = variablesToFFT(:,iVariable)
      jVariable = jVariable + 1
   end do
   write(*,*) "nVariablesToFFT=", grid%nVariablesToFFT, " variablesToFFT: ", grid%variablesToFFT

   do iStructureFactor = 1, nFactors
   do iPair = 1, 2
      do iVariable = 1, nVariablesToFFT

         if ( all(structureFactorPairs(:, iPair, iStructureFactor) == grid%variablesToFFT(:, iVariable)) ) &
            grid%structureFactorPairs(iPair, iStructureFactor) = iVariable

      end do
   end do
   end do
   write(*,*) "grid%structureFactorPairs=", grid%structureFactorPairs

   ! Initialize FFTW internals
   ! If the special direction is not periodic, we pretend like we have done an FFT
   ! but instead of wavenumber index it is really cell index along the special direction!
   !-------------------------------------
   iTemp=grid%nCells(grid%axisToPrint)
   
   if(.not.grid%periodic) then
      if(grid%axisToPrint/=2) then
         stop "If periodic=F only axisToPrint=2 is supported"
      end if
      grid%nCells(2)=1 ! Do not do FFTs along the special direction
   end if   
   allocate(grid%staticFFTarray(grid%nCells(1), grid%nCells(2), grid%nCells(3)))      
   call FFTW_PlanDFT(grid%staticFFTplan, grid%nCells(1), grid%nCells(2), grid%nCells(3), &
      grid%staticFFTarray, grid%staticFFTarray, FFTW_FORWARD, FFTW_ESTIMATE)
   grid%nCells(grid%axisToPrint) = iTemp   

   ! Some book-keeping
   grid%minWavenumber = ceiling((-grid%nCells + 1) / 2.0_wp)
   grid%maxWavenumber = floor (grid%nCells / 2.0_wp)
   print *, "Min k = ", grid%minWavenumber
   print *, "Max k = ", grid%maxWavenumber   

end subroutine createStaticFactors

subroutine createDynamicFactors(grid, wavenumbersString)
   type (HydroGrid), target :: grid
   character (*), intent(in) :: wavenumbersString

   integer :: iDimension, iWavenumber

   ! For the dynamic factors, we select which selectedWavenumbers to record:
   if(grid%nWavenumbers <= 0) return

   allocate(grid%selectedWavenumbers(3, grid%nWavenumbers))
   allocate(grid%selectedWaveindices(3, grid%nWavenumbers))
   if(grid%nCells(3)==1) then ! We read only two indices per wavenumber here
      read(wavenumbersString, *) grid%selectedWavenumbers(1:2,1:grid%nWavenumbers)
      grid%selectedWavenumbers(3,1:grid%nWavenumbers)=0      
   else
      read(wavenumbersString, *) grid%selectedWavenumbers
   end if
   do iWavenumber = 1, grid%nWavenumbers
      do iDimension = 1, 3
      if( (iDimension/=2) .or. (grid%periodic)) then
         if ( (grid%selectedWavenumbers(iDimension, iWavenumber) < grid%minWavenumber(iDimension)) .or. &
              (grid%selectedWavenumbers(iDimension, iWavenumber) > grid%maxWavenumber(iDimension)) ) then
              write(*,*) "Error: Selected wavenumber = ", iDimension, iWavenumber , " out of possible range"
              stop
         end if
         grid%selectedWaveindices(iDimension, iWavenumber) = frequencyToArrayIndex( &
            grid%selectedWavenumbers(iDimension, iWavenumber), grid%nCells(iDimension) )
      else 
         grid%selectedWaveindices(iDimension, iWavenumber) = grid%selectedWavenumbers(iDimension, iWavenumber)
      end if
      end do
   end do

   ! And for each of the selected wavenumbers, we keep a history (record) for some number of steps,
   ! and then do a temporal FFT and average the result over time:
   allocate(grid%savedStructureFactors(grid%nSavedSnapshots, grid%nWavenumbers, grid%nVariablesToFFT))
   allocate(grid%meanDynamicFactors(grid%nSavedSnapshots, grid%nWavenumbers, grid%nVariablesToFFT))
   allocate(grid%dynamicFactors(grid%nSavedSnapshots, grid%nWavenumbers, grid%nStructureFactors))

   ! Prepare FFTWs internal stuff:
   allocate(grid%dynamicFFTarray(grid%nSavedSnapshots))      
   call FFTW_PlanDFT(grid%dynamicFFTplan, grid%nSavedSnapshots, &
      grid%dynamicFFTarray, grid%dynamicFFTarray, FFTW_FORWARD, FFTW_ESTIMATE)

   ! Some book-keeping
   grid%minWavefrequency = ceiling((-grid%nSavedSnapshots + 1) / 2.0_wp)
   grid%maxWavefrequency = floor (grid%nSavedSnapshots / 2.0_wp)

end subroutine

!------------------------------------------------------------------------
! Destruction
!------------------------------------------------------------------------
subroutine destroyHydroAnalysis (grid)
   type (HydroGrid), target :: grid
   ! UNFINISHED
   write(*,*) "WARNING: Destruction of HydroGrid not really implemented"
end subroutine destroyHydroAnalysis

!------------------------------------------------------------------------
! Reset/initialize various means:
!------------------------------------------------------------------------

subroutine resetHydroAnalysis(grid)
   type (HydroGrid), target :: grid

   grid%meanPrimitive = 0
   if(grid%storeConserved) grid%meanConserved = 0
   grid%variances = 0
   grid%iStep = 0
   grid%iSample = 0
   grid%iTimeSeries = 0

   if(grid%nStructureFactors>0) then 
      grid%structureFactors = 0.0_wp
      if(allocated(grid%meanStaticTransforms)) grid%meanStaticTransforms = 0.0_wp
      if(grid%estimateCovariances) grid%structureFactorsCovariance = 0.0_wp
   end if

   if(grid%nWavenumbers>0) then
      grid%dynamicFactors = 0.0_wp
      grid%meanDynamicFactors = 0.0_wp
   end if

end subroutine resetHydroAnalysis

!------------------------------------------------------------------------
! Update given a snapshot of the hydrodynamic field
!------------------------------------------------------------------------

! This update routine accepts conserved variables as input, i.e., rho1, rho2, jx, jy, etc.   
subroutine updateHydroAnalysisConserved (grid, density, current, energy, scalars)
   type (HydroGrid), target :: grid
   real (wp), intent(in) :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0 : grid%nSpecies - 1)
   real (wp), intent(in) :: current(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nDimensions, 0:grid%nFluids) 
   real (wp), intent(in), optional :: energy(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in), optional :: scalars(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nPassiveScalars)

   ! Local variables
   real (wp) :: cellMass(0:grid%nSpecies - 1)
   real (wp) :: cellCurrent(grid%nDimensions, 0:grid%nFluids), cellEnergy(0:grid%nFluids)
   real (wp), dimension(0:grid%nFluids) :: temperature
   real (wp) :: densityOfSpeciesN, totalHeatCapacity
   integer :: iDimension, i, j, k, iSpecies

   integer :: mIdx, jIdx1, jIdx2, eIdx, cIdx1, cIdx2, sIdx1, sIdx2

   mIdx = grid%mIdx
   jIdx1 = grid%jIdx1
   jIdx2 = grid%jIdx2
   eIdx = grid%eIdx
   cIdx1 = grid%cIdx1
   cIdx2 = grid%cIdx2
   sIdx1 = grid%sIdx1
   sIdx2 = grid%sIdx2

   grid%iStep = grid%iStep + 1

   ! We need to convert all of the conserved to primitive variables in each cell
   do k = 1, grid%nCells(3)
   do j = 1, grid%nCells(2)
   do i = 1, grid%nCells(1)

      !grid%primitive(i, j, k, :, :) = 0 ! For debugging only

      ! Since we need to divide by density we should guard against zero density here
      grid%primitive(i, j, k, mIdx, 0:grid%nFluids) = MAX(8*TINY(density), density(i, j, k, 0:grid%nFluids))

      ! Velocities:
      do iDimension = 1, grid%nDimensions
         grid%primitive(i, j, k, jIdx1 + iDimension - 1, 0:grid%nFluids) = &
            current(i, j, k, iDimension, 0:grid%nFluids) / grid%primitive(i, j, k, mIdx, 0:grid%nFluids)
      end do

      ! Temperatures:
      cellMass(0) = grid%primitive(i, j, k, mIdx, 0) ! This is guaranteed to be strictly positive
      cellMass(1:grid%nSpecies-1) = density(i, j, k, 1:grid%nSpecies-1)
      cellCurrent = grid%primitive(i, j, k, jIdx1:jIdx2, 0:grid%nFluids)
      if(present(energy)) then
         cellEnergy = energy(i, j, k, 0:grid%nFluids)
         call energyToTemperature(grid, cellMass, cellCurrent, cellEnergy, temperature)
      else
         temperature = 1.0_wp
      end if   
      grid%primitive(i, j, k, eIdx, :) = temperature

      ! Concentrations:
      if (grid%isSingleFluid) then
         grid%primitive(i, j, k, cIdx1:cIdx2, 0) = density(i, j, k, 1:(cIdx2-cIdx1+1)) / grid%primitive(i, j, k, mIdx, 0)
      else if(cIdx2>=cIdx1) then ! We are storing concentrations as the last primitive variable
         grid%primitive(i, j, k, cIdx1, 1:grid%nFluids) = density(i, j, k, 1:grid%nFluids) / grid%primitive(i, j, k, mIdx, 0)
         grid%primitive(i, j, k, cIdx1, 0) = 1.0_wp - sum(grid%primitive(i, j, k, cIdx1, 1:grid%nFluids))
      end if

      ! Passive scalars (direct copy)
      if(present(scalars)) then
         grid%primitive(i, j, k, sIdx1:sIdx2, 0) = scalars(i, j, k, 1:(sIdx2-sIdx1+1))
      else
         grid%primitive(i, j, k, sIdx1:sIdx2, 0) = 0.0_wp   
      end if  

      if(grid%storeConserved) then
         ! We need the mean conserved variables to calculate unbiased means of the primitive variables
         grid%meanConserved(i, j, k, mIdx, :)         = ( (grid%iStep - 1) * grid%meanConserved(i, j, k, mIdx, :)        + &
            density(i, j, k, 0:grid%nFluids) ) / grid%iStep
         grid%meanConserved(i, j, k, jIdx1:jIdx2, :)  = ( (grid%iStep - 1) * grid%meanConserved(i, j, k, jIdx1:jIdx2, :) + &
            current(i, j, k, :, :) ) / grid%iStep
         if(present(energy)) then   
            grid%meanConserved(i, j, k, eIdx, :)         = ( (grid%iStep - 1) * grid%meanConserved(i, j, k, eIdx, :)        + &
               energy(i, j, k, :) ) / grid%iStep
         end if      
         if (grid%isSingleFluid) then
            ! The conserved variables include the partial densities
            grid%meanConserved(i, j, k, cIdx1:cIdx2, 0) = ( (grid%iStep-1) * grid%meanConserved(i, j, k, cIdx1:cIdx2, 0) + &
               density(i, j, k, 1:(cIdx2-cIdx1+1)) ) / grid%iStep
         end if
      end if   

   end do
   end do
   end do

   call updateHydroAnalysis (grid)      

end subroutine updateHydroAnalysisConserved

! This update routine accepts primitive variables as input, i.e., rho1, rho2, jx, jy, etc.   
subroutine updateHydroAnalysisPrimitive (grid, velocity, density, temperature, concentration, scalars, velocity_strided)
   type (HydroGrid), target :: grid
   ! Allow for alternative form of the storage format to avoid copy in/out for external codes
   real (wp), intent(in), optional :: &
      velocity(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nDimensions, 0:grid%nFluids), &
      velocity_strided(grid%nDimensions, grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)      
   real (wp), intent(in), optional :: density(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)   
   real (wp), intent(in), optional :: temperature(grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids)
   real (wp), intent(in), optional :: concentration(grid%nCells(1), grid%nCells(2), grid%nCells(3), 1:grid%nSpecies-1)
      ! Concentration is needed if treated as single fluid, otherwise it will be calculated
   real (wp), intent(in), optional :: scalars(grid%nCells(1), grid%nCells(2), grid%nCells(3), grid%nPassiveScalars)

   ! Local variables
   real (wp) :: cellMass(0:grid%nFluids)
   real (wp) :: densityOfSpeciesN
   integer :: iDimension, i, j, k, iSpecies

   integer :: mIdx, jIdx1, jIdx2, eIdx, cIdx1, cIdx2, sIdx1, sIdx2

   mIdx = grid%mIdx
   jIdx1 = grid%jIdx1
   jIdx2 = grid%jIdx2
   eIdx = grid%eIdx
   cIdx1 = grid%cIdx1
   cIdx2 = grid%cIdx2
   sIdx1 = grid%sIdx1
   sIdx2 = grid%sIdx2

   grid%iStep = grid%iStep + 1
   if(grid%storeConserved) then
      write(0,*) "WARNING: It is not consisent to have storeConserved=T when using primitive variables!"
      write(0,*) "WARNING: Switching grid%storeConserved=F"
      grid%storeConserved=.false.
   end if   

   !if(present(scalars)) write(*,*) sIdx1,sIdx2, " Scalars=", real(scalars(:, :, :, 1:(sIdx2-sIdx1+1)))
   if(.false.) write(*,*) "HydroGrid present=", present(velocity), present(density), &
      present(concentration), present(temperature), present(scalars)

   ! We need to convert all of the conserved to primitive variables in each cell
   do k = 1, grid%nCells(3)
   do j = 1, grid%nCells(2)
   do i = 1, grid%nCells(1)

      if(present(velocity)) then
         do iDimension = 1, grid%nDimensions
            grid%primitive(i, j, k, jIdx1 + iDimension - 1, 0:grid%nFluids) = velocity(i, j, k, iDimension, 0:grid%nFluids)
         end do
      else if(present(velocity_strided)) then
         do iDimension = 1, grid%nDimensions
            grid%primitive(i, j, k, jIdx1 + iDimension - 1, 0:grid%nFluids) = velocity_strided(iDimension, i, j, k, 0:grid%nFluids)
         end do
      else
         do iDimension = 1, grid%nDimensions
            grid%primitive(i, j, k, jIdx1 + iDimension - 1, 0:grid%nFluids) = 0.0_wp ! Default value
         end do
      end if
     
      if(present(density)) then
         grid%primitive(i, j, k, mIdx, 0:grid%nFluids) = density(i, j, k, 0:grid%nFluids)
      else
         grid%primitive(i, j, k, mIdx, 0:grid%nFluids) = 1.0_wp ! Default value   
      end if   

      if(present(temperature)) then
         grid%primitive(i, j, k, eIdx, 0:grid%nFluids) = temperature(i, j, k, 0:grid%nFluids)
      else
         grid%primitive(i, j, k, eIdx, 0:grid%nFluids) = 1.0_wp ! Default value   
      end if   

      ! Concentrations:
      if (grid%isSingleFluid) then
         if(.not.present(concentration)) then
            if(grid%nSpecies>1) then
               stop "Concentration must be present in updateHydroAnalysisPrimitive ifSingleFluid for mixtures!"
            end if   
         else   
            grid%primitive(i, j, k, cIdx1:cIdx2, 0) = concentration(i, j, k, 1:grid%nSpecies-1)
         end if   
      else if(cIdx2>=cIdx1) then ! We are storing concentrations as the last primitive variable
         grid%primitive(i, j, k, cIdx1, 1:grid%nFluids) = density(i, j, k, 1:grid%nFluids) / density(i, j, k, 0)        
         grid%primitive(i, j, k, cIdx1, 0) = 1.0_wp - sum(grid%primitive(i, j, k, cIdx1, 1:grid%nFluids))
      end if

      ! Passive scalars (direct copy)
      if(present(scalars)) then
         grid%primitive(i, j, k, sIdx1:sIdx2, 0) = scalars(i, j, k, 1:(sIdx2-sIdx1+1))
      else
         grid%primitive(i, j, k, sIdx1:sIdx2, 0) = 0.0_wp   
      end if  

   end do
   end do
   end do

   call updateHydroAnalysis (grid)
   
end subroutine updateHydroAnalysisPrimitive

! This update routine is meant to be used in staggered codes, even though our analysis kind of assumes cell-centered variables
subroutine updateHydroAnalysisStaggered (grid, nGhost, vx, vy, vz, nGhostScalar, density)
   type (HydroGrid), target :: grid
   integer, intent(in) :: nGhost
   ! Sometimes the velocities may be stored separately, on different grids, as in staggered codes:
   real (wp), intent(in), dimension(grid%nCells(1)+nGhost, grid%nCells(2), grid%nCells(3), 0:grid%nFluids) :: vx
   real (wp), intent(in), dimension(grid%nCells(1), grid%nCells(2)+nGhost, grid%nCells(3), 0:grid%nFluids) :: vy
   real (wp), intent(in), dimension(grid%nCells(1), grid%nCells(2), grid%nCells(3)+nGhost, 0:grid%nFluids), optional :: vz
   !real (wp), intent(in), dimension(1-nGhost:grid%nCells(1), grid%nCells(2), grid%nCells(3), 0:grid%nFluids) :: vx
   !real (wp), intent(in), dimension(grid%nCells(1), 1-nGhost:grid%nCells(2), grid%nCells(3), 0:grid%nFluids) :: vy
   !real (wp), intent(in), dimension(grid%nCells(1), grid%nCells(2), 1-nGhost:grid%nCells(3), 0:grid%nFluids), optional :: vz   
   integer, intent(in) :: nGhostScalar ! Cannot be optional here as used in specification expression below
   real (wp), intent(in), optional :: density(1-nGhostScalar:grid%nCells(1)+nGhostScalar, &
                                              1-nGhostScalar:grid%nCells(2)+nGhostScalar, &
                  1-merge(nGhostScalar,0,grid%nDimensions>2):grid%nCells(3)+merge(nGhostScalar,0,grid%nDimensions>2), &
                                              0:grid%nFluids)

   ! Local variables
   real (wp) :: cellMass(0:grid%nFluids)
   real (wp) :: densityOfSpeciesN
   integer :: iDimension, i, j, k, iSpecies

   integer :: mIdx, jIdx1, jIdx2, eIdx, cIdx1, cIdx2, sIdx1, sIdx2
   
   !write(*,*) "TEST:", "nCells=",grid%nCells, " nGhost=", nGhost

   mIdx = grid%mIdx
   jIdx1 = grid%jIdx1
   jIdx2 = grid%jIdx2
   eIdx = grid%eIdx
   cIdx1 = grid%cIdx1
   cIdx2 = grid%cIdx2
   sIdx1 = grid%sIdx1
   sIdx2 = grid%sIdx2

   grid%iStep = grid%iStep + 1
   if(grid%storeConserved) then
      write(0,*) "WARNING: It is not consisent to have storeConserved=T when using primitive variables!"
      write(0,*) "WARNING: Switching grid%storeConserved=F"
      grid%storeConserved=.false.
   end if   

   !if(present(density)) write(*,*) sIdx1,sIdx2, " Scalars=", real(density(:, :, :, 1:(sIdx2-sIdx1+1)))

   ! We need to convert all of the conserved to primitive variables in each cell
   do k = 1, grid%nCells(3)
   do j = 1, grid%nCells(2)
   do i = 1, grid%nCells(1)

      grid%primitive(i, j, k, jIdx1 - 1 + 1, 0:grid%nFluids) = vx(i, j, k, 0:grid%nFluids)
      grid%primitive(i, j, k, jIdx1 - 1 + 2, 0:grid%nFluids) = vy(i, j, k, 0:grid%nFluids)
      if(grid%nDimensions>2) then ! vz must be present if 3D         
         grid%primitive(i, j, k, jIdx1 - 1 + 3, 0:grid%nFluids) = vz(i, j, k, 0:grid%nFluids)
      end if   
     
      if(present(density)) then
         grid%primitive(i, j, k, mIdx, 0:grid%nFluids) = density(i, j, k, 0:grid%nFluids)
      else
         grid%primitive(i, j, k, mIdx, 0:grid%nFluids) = 1.0_wp ! Default value   
      end if   

      ! Default values for the rest:
      grid%primitive(i, j, k, eIdx, 0:grid%nFluids) = 1.0_wp ! Default value
      grid%primitive(i, j, k, cIdx1:cIdx2, 0) = 0.0_wp
      grid%primitive(i, j, k, sIdx1:sIdx2, 0) = 0.0_wp   

   end do
   end do
   end do

   call updateHydroAnalysis (grid)
   
end subroutine updateHydroAnalysisStaggered

subroutine updateHydroAnalysis (grid)
   type (HydroGrid), target :: grid
   
   real(wp) :: conc
   integer :: iDimension, i, j, k, iSpecies

   if((useRelativeVelocity>=0) .and. (grid%nFluids>1)) then
      do k = 1, grid%nCells(3)
      do j = 1, grid%nCells(2)
      do i = 1, grid%nCells(1)
         ! We need the mean concentration profile here:
         conc= grid%bottomConcentration + ((j-0.5_wp)/grid%nCells(2))*(grid%topConcentration-grid%bottomConcentration)
         if(useRelativeVelocity==0) then
            ! Calculate (1-c_bar)*v12, where v12=v1-v2 and store it in v2
            ! And replace rho with c_bar*rho
            grid%primitive(i, j, k, grid%mIdx, 0) = conc * grid%primitive(i, j, k, grid%mIdx, 0)
            grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, grid%nFluids) = (1.0_wp-conc)* &
               (grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, 1) - &
               grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, grid%nFluids))
         else
            ! Calculate (1-c_bar)*c_bar*v12, where v12=v1-v2 and store it in v2
            grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, grid%nFluids) = conc*(1.0_wp-conc)* &
               (grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, 1) - &
               grid%primitive(i, j, k, grid%jIdx1:grid%jIdx2, grid%nFluids))
         end if      
      end do
      end do
      end do
   
   else if((useRelativeVelocity>=0) .and. (grid%nFluids==0) .and. (grid%nSpecies>2)) then
      ! There is no relative velocity in single-fluid models, so we need to improvise
            
      do k = 1, grid%nCells(3)
      do j = 1, grid%nCells(2)
      do i = 1, grid%nCells(1)
         ! We need the mean concentration profile here:
         conc= grid%bottomConcentration + ((j-0.5_wp)/grid%nCells(2))*(grid%topConcentration-grid%bottomConcentration)
         ! Here I replace c2 with drho1 / [c*(1-c)], where rho1 is calculated as rho*c
         grid%primitive(i, j, k, grid%cIdx2, grid%nFluids) = &
            ( grid%primitive(i, j, k, grid%mIdx, grid%nFluids) * grid%primitive(i, j, k, grid%cIdx1, grid%nFluids) ) / &
            (conc*(1.0_wp-conc))
      end do
      end do
      end do      
      
   end if

   call updateMeans(grid)
   call updateVariances(grid)
   call updateStructureFactors(grid)      

end subroutine


subroutine updateMeans(grid) ! Update the means of primitive variables
   type (HydroGrid), target :: grid
   
   ! We keep a running average of the mean that is correct at any point in time:
   grid%meanPrimitive = ( (grid%iStep - 1) * grid%meanPrimitive + grid%primitive ) / grid%iStep

end subroutine


subroutine updateVariances(grid)
   type (HydroGrid), target :: grid
   
   integer :: iVariance, variable1, species1, variable2, species2
   integer, dimension(nMaxDims) :: cell
   
   cell=grid%correlationCell

   do iVariance = 1, grid%nVariances
      species1  = grid%variancePairs(1, 1, iVariance)
      variable1 = abs(grid%variancePairs(2, 1, iVariance))
      species2  = grid%variancePairs(1, 2, iVariance)
      variable2 = abs(grid%variancePairs(2, 2, iVariance))
      
      !write(*,*) iVariance, variable1, species1, variable2, species2
      !write(*,*) "X=", grid%primitive(:, :, :, variable1, species1)

      ! Keep track of <AB>, for selected pairs of variables A and B
      if(grid%variancePairs(2, 1, iVariance)>0) then      
         grid%variances(:, :, :, iVariance) = ( (grid%iStep - 1) * grid%variances(:, :, :, iVariance) + &
            grid%primitive(:, :, :, variable1, species1) * grid%primitive(:, :, :, variable2, species2) ) / grid%iStep
      else
         grid%variances(:, :, :, iVariance) = ( (grid%iStep - 1) * grid%variances(:, :, :, iVariance) + &
            grid%primitive(cell(1), cell(2), cell(3), variable1, species1) * &
            grid%primitive(:, :, :, variable2, species2) ) / grid%iStep      
      end if      

   end do

end subroutine updateVariances   


subroutine updateStructureFactors(grid)
   type (HydroGrid), target :: grid

   integer :: iVariable, species, variable, iTemp
   integer :: iStructureFactor, jStructureFactor, variable1, variable2, vars(4)
   integer :: iWavenumber, kx, ky, kz, iCross, i, j, k, vel_index

   ! These will be used to index the structure factors, starting at the location of the negative wavenumbers.
   ! For example, if nCells(1) = 6
   ! the corresponding 1-D slice has components 0 1 2 3 -2 -1.
   ! Thus, iCell will start with the value 5 since that is the location of the 
   ! first negative frequency, which is -2 * (2*pi/L).  Both even and odd nCell values should work
   integer :: iCell, jCell, kCell
   
   integer, dimension(nMaxDims) :: mink, maxk, ijk
   real(wp), dimension(nMaxDims) :: kdrh
   character (nMaxCharacters) :: filenameBase=""

   !---------------------------------------
   ! Static structure factors
   if (grid%nStructureFactors <= 0) return

   mink=grid%minWavenumber
   maxk=grid%maxWavenumber

   do iVariable = 1, grid%nVariablesToFFT
      species = grid%variablesToFFT(1, iVariable)
      variable = grid%variablesToFFT(2, iVariable)
   
      if(.not.grid%periodic) then
         ! Do a series of two-dimensional FFTs along the special axis
         jCell = frequencyToArrayIndex(mink(2), grid%nCells(2))
         do j = mink(2), maxk(2)
            grid%staticFFTarray(:, 1, :) = grid%primitive(:, j-mink(2)+1, :, variable, species)
            call FFTW_Execute(grid%staticFFTplan)       
            grid%staticFourierTransforms(:, jCell, :, iVariable) = grid%staticFFTarray(:, 1, :)
            jCell = jCell + 1 ; if (jCell > grid%nCells(2)) jCell = 1
         end do
      else ! Do one big 3D FFT along all axes
         grid%staticFFTarray = grid%primitive(:, :, :, variable, species)      
         call FFTW_Execute(grid%staticFFTplan)       
         grid%staticFourierTransforms(:, :, :, iVariable) = grid%staticFFTarray
      end if   

      ! We may need to correct for staggering of the velocities by multiplying by the phase shifts exp(-I*k*dx/2)
      vel_index = variable - grid%jIdx1 + 1
      if( (abs(grid%staggeredVelocities)>0) .and. (vel_index>=1) .and. (vel_index<=grid%nDimensions)) then      
         !write(*,*) "Correcting staggering of variable ", variable, " index=", vel_index
         
         iCell = frequencyToArrayIndex(mink(1), grid%nCells(1))
         jCell = frequencyToArrayIndex(mink(2), grid%nCells(2))
         kCell = frequencyToArrayIndex(mink(3), grid%nCells(3))
         do k = mink(3), maxk(3)
         do j = mink(2), maxk(2)
         do i = mink(1), maxk(1)
            ijk = (/i,j,k/)
            
            kdrh = grid%staggeredVelocities * ijk * pi / grid%nCells ! The sign could be +1 or -1
            if(.not.grid%periodic) kdrh(2)=0.0_wp ! No phase shift along y
            
            grid%staticFourierTransforms(iCell, jCell, kCell, iVariable) = exp(cmplx(0.0_wp, kdrh(vel_index), wp)) * &
            grid%staticFourierTransforms(iCell, jCell, kCell, iVariable)               
               
         iCell = iCell + 1 ; if (iCell > grid%nCells(1)) iCell = 1
         end do
         jCell = jCell + 1 ; if (jCell > grid%nCells(2)) jCell = 1
         end do
         kCell = kCell + 1 ; if (kCell > grid%nCells(3)) kCell = 1
         end do
      
      end if

      ! Keep a runnning average of the static transforms:
      if(allocated(grid%meanStaticTransforms)) grid%meanStaticTransforms(:, :, :, iVariable) = &
         ( (grid%iStep - 1) * grid%meanStaticTransforms(:, :, :, iVariable) + &
         grid%staticFourierTransforms(:, :, :, iVariable) ) / grid%iStep

   end do

   do iStructureFactor = 1, grid%nStructureFactors
      variable1 = grid%structureFactorPairs(1, iStructureFactor)
      variable2 = grid%structureFactorPairs(2, iStructureFactor)
      ! Keep track of <A_hat*conj(B_hat)> for selected pairs of variables A and B:
      grid%structureFactors(:, :, :, iStructureFactor) = &
         ( (grid%iStep - 1) * grid%structureFactors(:, :, :, iStructureFactor) + &
         grid%staticFourierTransforms(:, :, :, variable1) * conjg(grid%staticFourierTransforms(:, :, :, variable2)) &
         ) / grid%iStep
   end do

   if(grid%estimateCovariances) then ! Correlations between different structure factors
      ! This might be useful in variance reduction  
      ! Here we only look at the real part of the structure factors!
      
      iCross=0
      do iStructureFactor = 1, grid%nStructureFactors
      do jStructureFactor = 1, grid%nStructureFactors
      if(iStructureFactor<=jStructureFactor) then
         iCross=iCross+1
         vars(1) = grid%structureFactorPairs(1, iStructureFactor)
         vars(2) = grid%structureFactorPairs(2, iStructureFactor)
         vars(3) = grid%structureFactorPairs(1, jStructureFactor)
         vars(4) = grid%structureFactorPairs(2, jStructureFactor)
         grid%structureFactorsCovariance(:, :, :, iCross) = &
            ( (grid%iStep - 1) * grid%structureFactorsCovariance(:, :, :, iCross) + &
              real(grid%staticFourierTransforms(:, :, :, vars(1)) * conjg(grid%staticFourierTransforms(:, :, :, vars(2)))) * &
              real(grid%staticFourierTransforms(:, :, :, vars(3)) * conjg(grid%staticFourierTransforms(:, :, :, vars(4)))) &
            ) / grid%iStep
      end if
      end do
      end do      
   end if   

   !---------------------------------------
   ! Dynamic structure factors
   if (grid%nWavenumbers <= 0) return

   if (grid%iSample >= abs(grid%nSavedSnapshots)) then ! Memory of past history is full, do temporal FFT now
      grid%iSample = 0
      grid%iTimeSeries = grid%iTimeSeries + 1

      do iVariable = 1, grid%nVariablesToFFT         
         do iWavenumber = 1, grid%nWavenumbers
            grid%dynamicFFTarray = grid%savedStructureFactors(:, iWavenumber, iVariable)
            call FFTW_Execute(grid%dynamicFFTplan)
            grid%savedStructureFactors(:, iWavenumber, iVariable) = grid%dynamicFFTarray ! Use this as temporary storage

            ! Also keep track of the mean of the FT:
            grid%meanDynamicFactors(:, iWavenumber, iVariable) = &
               ( (grid%iTimeSeries - 1) * grid%meanDynamicFactors(:, iWavenumber, iVariable) + &
               grid%savedStructureFactors(:, iWavenumber, iVariable) ) / grid%iTimeSeries
         end do
      end do

      ! Keep a running time average of the dynamic structure factors
      do iStructureFactor = 1, grid%nStructureFactors
         variable1 = grid%structureFactorPairs(1, iStructureFactor)
         variable2 = grid%structureFactorPairs(2, iStructureFactor)

         ! Recall shape is grid%dynamicFactors(grid%nSavedSnapshots, grid%nWavenumbers, grid%nStructureFactors)
         if(grid%correlationWavenumber>0) then
            do iWavenumber = 1, grid%nWavenumbers
               grid%dynamicFactors(:, iWavenumber, iStructureFactor) = &
               ( (grid%iTimeSeries - 1) * grid%dynamicFactors(:, iWavenumber, iStructureFactor) + &
               grid%savedStructureFactors(:, grid%correlationWavenumber, variable1) * &
               conjg(grid%savedStructureFactors(:, iWavenumber, variable2)) ) / grid%iTimeSeries
            end do          
         else
            grid%dynamicFactors(:, :, iStructureFactor) = &
               ( (grid%iTimeSeries - 1) * grid%dynamicFactors(:, :, iStructureFactor) + &
               grid%savedStructureFactors(:, :, variable1) * conjg(grid%savedStructureFactors(:, :, variable2)) &
               ) / grid%iTimeSeries
         end if      
      end do

   end if

   ! Save the current snapshot in the record (history):
   grid%iSample = grid%iSample + 1
   do iVariable = 1, grid%nVariablesToFFT
      do iWavenumber = 1, grid%nWavenumbers
         kx = grid%selectedWaveindices(1, iWavenumber)
         ky = grid%selectedWaveindices(2, iWavenumber)
         kz = grid%selectedWaveindices(3, iWavenumber)
         grid%savedStructureFactors(grid%iSample, iWavenumber, iVariable) = grid%staticFourierTransforms(kx, ky, kz, iVariable)
      end do
   end do

   if(.false.) then ! Temporary: Save a history of S(k)
      filenameBase = trim(grid%outputFolder) // "/" // trim(grid%filePrefix) // ".S_k.history.dat"
      write(*,*) "Writing S(k) history to file ", trim(filenameBase)
      open (551, file = trim(filenameBase), status = "unknown", action = "write", position="append")

      write(551,'(1000g17.9)') &
         abs  (grid%savedStructureFactors(grid%iSample, :, :)), &
         real (grid%savedStructureFactors(grid%iSample, :, :)), &
         aimag(grid%savedStructureFactors(grid%iSample, :, :))
         
      close(551)
   end if   

end subroutine updateStructureFactors


!------------------------------------------------------------------------
! File I/O
!------------------------------------------------------------------------      

subroutine writeToFiles(grid, id)
   ! For now this only writes one-dimensional profiles along a given dimension (axisToPrint)
   ! This saves the means, variances, and structure factors
   type (HydroGrid), target :: grid
   integer, optional :: id ! An additional integer to append to file names

   integer :: iVariance, species1, species2, variable1, variable2
   integer, dimension(nMaxDims) :: cell
   
   character (nMaxCharacters) :: filenameBase=""
   character(25) :: id_string
   
   cell=grid%correlationCell

   ! First, subtract the means from the variances to get covariances:
   ! If we only have one sample we cannot evaluate the means so we should not subtract them
   if(grid%iStep>1) then
   do iVariance = 1, grid%nVariances
      species1  = grid%variancePairs(1, 1, iVariance)
      variable1 = abs(grid%variancePairs(2, 1, iVariance))
      species2  = grid%variancePairs(1, 2, iVariance)
      variable2 = abs(grid%variancePairs(2, 2, iVariance))

      !write(*,*) "ID=", iVariance, species1, species2, variable1, variable2
      !write(*,*) " var=", grid%variances(:,:,:,iVariance)
      !write(*,*) "mean=", grid%meanPrimitive(:,:,:,variable1,species1)
      
      if(grid%variancePairs(2, 1, iVariance)>0) then  
         grid%variances(:,:,:,iVariance) = grid%variances(:,:,:,iVariance) - &
           grid%meanPrimitive(:,:,:,variable1,species1) * grid%meanPrimitive(:,:,:,variable2,species2)
      else    
         grid%variances(:,:,:,iVariance) = grid%variances(:,:,:,iVariance) - &
           grid%meanPrimitive(cell(1),cell(2),cell(3),variable1,species1) * &
           grid%meanPrimitive(:,:,:,variable2,species2)
         if((species1==species2).and.(variable1==variable2)) then ! Normalize to eliminate scaling
            ! Normalize by the variance to get correlation coefficient:
            grid%variances(:,:,:,iVariance) = grid%variances(:,:,:,iVariance) / &
               grid%variances(cell(1),cell(2),cell(3),iVariance) / grid%structFactMultiplier
            ! Eliminate the delta function peak at the origin:   
            grid%variances(cell(1),cell(2),cell(3),iVariance) = 0
         end if   
           
      end if
              
   end do
   end if

   filenameBase = trim(grid%outputFolder) // "/" // trim(grid%filePrefix)

   if(present(id)) then
      write(id_string,"(I8.8)") id
      filenameBase = trim(filenameBase) // "." // trim(ADJUSTL(id_string))
   end if   

   write(*,*) "Writing HydroGrid files for  ", trim(filenameBase), &
      " from # of samples / timeseries =", grid%iStep, grid%iTimeSeries

   call writeMeansAndVariances(grid, trim(filenameBase))
   call writeStructureFactors(grid, trim(filenameBase))

   call resetHydroAnalysis(grid) ! Some of the arrays were changed before writing

end subroutine

subroutine writeMeansAndVariances(grid, filenameBase)
   ! For now this only writes one-dimensional profiles along a given dimension (axisToPrint)
   ! This saves the means, variances, and structure factors
   type (HydroGrid), target :: grid
   character(len=*) :: filenameBase

   integer, parameter :: unitOffset = 330
   integer :: meanFile = unitOffset+1, instantMeanFile = unitOffset+2, varianceFile = unitOffset+3

   integer :: otherAxes(2)
   real(wp), dimension(grid%nCells(grid%axisToPrint), grid%nVariables, 0:grid%nFluids)  :: instantMeans
   real(wp), dimension(grid%nCells(grid%axisToPrint), grid%nVariances) :: variances

   ! Mean conserved variables:
   real(wp), dimension(grid%nCells(grid%axisToPrint), 0:grid%nSpecies-1)  :: meanDensity
      ! We always include the partial densities as well
   real(wp), dimension(grid%nCells(grid%axisToPrint), 0:grid%nFluids)  :: meanTemperature, meanEnergy
   ! Unbiased primitive means:
   real(wp), dimension(grid%nCells(grid%axisToPrint), grid%nDimensions, 0:grid%nFluids)  :: meanVelocity

   real(wp) :: meanConcentration(0:grid%nSpecies-1)
   integer :: i, j, k, iSpecies, species, species1, species2, iVariable, variable, variable1, variable2, iVariance
   integer :: iCell, jCell, kCell, iDimension

   integer :: mIdx, jIdx1, jIdx2, eIdx, cIdx1, cIdx2, sIdx1, sIdx2, axis
   mIdx = grid%mIdx
   jIdx1 = grid%jIdx1
   jIdx2 = grid%jIdx2
   eIdx = grid%eIdx
   cIdx1 = grid%cIdx1
   cIdx2 = grid%cIdx2
   sIdx1 = grid%sIdx1
   sIdx2 = grid%sIdx2      
   axis=grid%axisToPrint

   if((axis>3).or.(axis<1)) then
      stop "Must select a special axis for outputting statistics"
   end if
   otherAxes = (/2, 1/)
   where (otherAxes == axis) otherAxes = 3  
   if (otherAxes(1) < otherAxes(2)) otherAxes = otherAxes((/2, 1/))

   ! Spatial means of instantaneous values, e.g., <v>=<j/rho> (these have a 1/N bias)     
   !-----------------------     
   instantMeans = sum(sum(grid%meanPrimitive, otherAxes(1)), otherAxes(2)) / &
      ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) )
   ! One may also calculate instantaneous means using the spatially-averaged instantaneous values
   ! These will still have a bias, but smaller if spatial averaging over multiple cells is performed

   !-----------------------
   ! Now actually write the mean profiles to files
   write(*,*) "Writing instantaneous (biased) means to file ", trim(filenameBase) // ".means.inst.dat"
   open (instantMeanFile, file = trim(filenameBase) // ".means.inst.dat", status = "unknown", action = "write")

   do iSpecies=0, grid%nSpecies-1
      write(instantMeanFile, '(A,I0,100(g17.9))') "# species/mean=", iSpecies, &
         sum(instantMeans(1:grid%nCells(axis), :, 0), dim=1)/grid%nCells(axis)
   end do   
   
   if (grid%isSingleFluid) then      
      call writeHeaderSingleFluid(instantMeanFile)      
      do iCell = 1, grid%nCells(axis)
         write(instantMeanFile, '(100(g17.9))') (iCell-0.5_wp)*grid%systemLength(axis)/grid%nCells(axis), &
            instantMeans(iCell, :, 0)
      end do
   else
      call writeHeaderMultiFluid(instantMeanFile)
      do iCell = 1, grid%nCells(axis)
         write(instantMeanFile, '(100(g17.9))')  (iCell-0.5_wp)*grid%systemLength(axis)/grid%nCells(axis),  &
            (instantMeans(iCell, :, iSpecies), iSpecies=0, grid%nSpecies-1)
      end do
   end if

   close (instantMeanFile)
      
   if(grid%writeSnapshotVTK) then
      call writeSnapshotToVTK()
   end if   
   
   if(grid%writeMeansVTK) then
      call writeVelocityToVTK()
   end if   

   ! Now unbiased means, e.g., <v>=<j>/<rho>, where the average is both temporal and spatial to ensure minimal fluctuations
   !-----------------------     
   if(grid%storeConserved) then
   
      meanDensity(:, 0:grid%nFluids)  = sum(sum(grid%meanConserved(:,:,:,mIdx,:), otherAxes(1)), otherAxes(2)) / &
         ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) ) ! <rho>
      if(grid%isSingleFluid) then
         meanDensity(:, 1:(cIdx2-cIdx1+1))  = sum(sum(grid%meanConserved(:,:,:,cIdx1:cIdx2,0), otherAxes(1)), otherAxes(2)) / &
            ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) ) ! <rho1>, <rho2>, etc.
      end if         

      ! <v>=<j>/<rho>
      do iDimension = 1, grid%nDimensions
         meanVelocity(:, iDimension, :) = &
            sum(sum(grid%meanConserved(:, :, :, jIdx1 + iDimension - 1, :), otherAxes(1)), otherAxes(2)) &
               / ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) ) / meanDensity(:, 0:grid%nFluids)
      end do

      meanEnergy = sum(sum(grid%meanConserved(:,:,:,eIdx,:), otherAxes(1)), otherAxes(2)) / &
         ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) ) ! <e>

      do iCell = 1, grid%nCells(axis) ! Convert to temperature <T>
         call energyToTemperature(grid, density=meanDensity(iCell,:), velocity=meanVelocity(iCell, :, :), &
            energy=meanEnergy(iCell, :), temperature=meanTemperature(iCell,:))
      end do

      write(*,*) "Writing unbiased means to file ", trim(filenameBase) // ".means.dat"
      open (meanFile,        file = trim(filenameBase) // ".means.dat",      status = "unknown", action = "write")

      if (grid%isSingleFluid) then
         call writeHeaderSingleFluid(meanFile)           
      else
         call writeHeaderMultiFluid(meanFile)
      end if         
      do iCell = 1, grid%nCells(axis)
         meanConcentration(1:grid%nSpecies-1) = meanDensity(iCell,1:grid%nSpecies-1) / meanDensity(iCell,0) ! Unbiased estimate
         meanConcentration(0) = 1.0_wp - sum(meanConcentration(1:grid%nSpecies-1))
         if (grid%isSingleFluid) then
            write(meanFile, '(100(g17.9))') (iCell-0.5_wp)*grid%systemLength(axis)/grid%nCells(axis),  &
               meanDensity(iCell,0), meanVelocity(iCell, :, 0), meanTemperature(iCell,0), meanConcentration(1:grid%nSpecies-1)
         else      
            write(meanFile, '(100(g17.9))') (iCell-0.5_wp)*grid%systemLength(axis)/grid%nCells(axis), &
               (meanDensity(iCell,iSpecies), meanVelocity(iCell, :, iSpecies), meanTemperature(iCell,iSpecies), &
                meanConcentration(iSpecies), iSpecies=0, grid%nSpecies-1)
         end if
      end do         

      close (meanFile)

   end if      

   if(grid%nVariances>0) then

      ! Account for a normalization of the noise:
      grid%variances = grid%variances * grid%structFactMultiplier
   
      !-----------------------
      ! Mean variances (always have a bias)
      variances = sum(sum(grid%variances, otherAxes(1)), otherAxes(2)) / ( grid%nCells(otherAxes(1)) * grid%nCells(otherAxes(2)) )
      write(*,*) "Writing averaged variances to file ", trim(filenameBase) // ".variances.dat"
      open (varianceFile, file = trim(filenameBase) // ".variances.dat",  &
         status = "unknown", action = "write")

      write(varianceFile, '(A,100(g17.9))') "# variances=", sum(variances(1:grid%nCells(axis), :), dim=1)/grid%nCells(axis)
      write(varianceFile, '(A,100(g17.9))') "# |variances|=", sum(abs(variances(1:grid%nCells(axis), :)), dim=1)/grid%nCells(axis)

      write(varianceFile, "(A, 100('[(',I0,I0,'),(',I0,I0,')]'))") "# position, variances of: ", &
         (grid%variancePairs(1:2, 1:2, iVariance), iVariance=1,grid%nVariances)
      write(varianceFile, '(100(A,I0))') "# nCells=", grid%nCells(axis), &
         " nVariances=", grid%nVariances, " nVariables=", grid%nVariables, &
         " nVelocities=", grid%nDimensions, " nConcentrations=", grid%nSpecies-1, " nScalars=", grid%nPassiveScalars

      do iCell = 1, grid%nCells(axis)
         write(varianceFile, '(100(g17.9))') (iCell-0.5_wp)*grid%systemLength(axis)/grid%nCells(axis), variances(iCell, :)
      end do
      close (varianceFile)

      if(grid%writeVariancesVTK) then
         call writeVariancesToVTK()
      end if
   end if
   
contains

   subroutine writeHeaderSingleFluid(file)
      integer, intent(in) :: file
      write(file, '(A)') "# position, rho, v, T, c, s"
      write(file, '(100(A,I0))') "# nCells=", grid%nCells(axis), " nVariables=", grid%nVariables, &
         " nVelocities=", grid%nDimensions, " nConcentrations=", grid%nSpecies-1, " nScalars=", grid%nPassiveScalars       
   end subroutine

   subroutine writeHeaderMultiFluid(file)
      integer, intent(in) :: file
      write(file, '(A)') "# position, rho, v, T, c(nSpecies) and then (rho,v,T,c) for species 1:nSpecies-1"
      write(file, '(100(A,I0))') "# nCells=", grid%nCells(axis), " nSpecies=", grid%nSpecies, &
         " nVariables=", grid%nVariables, " nVelocities=", grid%nDimensions     
   end subroutine

   subroutine writeSnapshotToVTK()
      ! Write a snapshot of the instantaneous fields
      integer :: mesh_dims(3), dim, iVariance
      character(len=16), dimension(max(6,grid%nVariables)), target :: varnames

      character(len=nMaxCharacters), target :: filename
      
      real(wp), dimension(:,:,:,:), allocatable, target :: velocity
      integer :: nvars

      filename = trim(filenameBase) // ".snapshot.vtk"
      write(*,*) "Writing instantaneous single-fluid variables to file ", trim(filename)
      filename = trim(filename) // C_NULL_CHAR

      mesh_dims=grid%nCells(1:3)
      if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)

      if(grid%nDimensions==0) then ! Don't plot velocity, plot scalars only
         nvars=grid%nVariables-grid%nDimensions
         varnames(1) = "Density" // C_NULL_CHAR
         varnames(2) = "Temperature" // C_NULL_CHAR
         varnames(3) = "Scalar1" // C_NULL_CHAR
         varnames(4) = "Scalar2" // C_NULL_CHAR
         varnames(5) = "Scalar3" // C_NULL_CHAR
         varnames(6:) = "Scalar" // C_NULL_CHAR
         CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
            ub=0_c_int, dims=mesh_dims+1, &
            x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
            y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
            z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
            nvars=nvars, vardim=(/(1, dim=1,nvars)/), &
            centering=(/(0, dim=1,nvars)/), &
            varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,nvars) /), &
            vars=(/ C_LOC(grid%primitive(1,1,1,grid%mIdx,0)), &
                   (C_LOC(grid%primitive(1,1,1,grid%eIdx+dim,0)), dim=0,grid%nVariables-grid%eIdx) /))
      else

         ! We always make the velocity be three dimensional to plot correctly in visit as a vector
         allocate(velocity(3, grid%nCells(1), grid%nCells(2), grid%nCells(3)))
         do dim=1, grid%nDimensions
            velocity(dim, :, :, :) = grid%primitive(:, :, :, grid%jIdx1 + dim - 1 , 0)
         end do   
         velocity(grid%nDimensions+1 : 3, :, :, :) = 0.0_wp

         nvars=grid%nVariables+1-grid%nDimensions
         varnames(1) = "Density" // C_NULL_CHAR
         varnames(2) = "Velocity" // C_NULL_CHAR
         varnames(3) = "Temperature" // C_NULL_CHAR
         varnames(4) = "Scalar1" // C_NULL_CHAR
         varnames(5) = "Scalar2" // C_NULL_CHAR
         varnames(6) = "Scalar3" // C_NULL_CHAR
         varnames(7:) = "Scalar" // C_NULL_CHAR
         CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
            ub=0_c_int, dims=mesh_dims+1, &
            x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
            y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
            z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
            nvars=nvars, vardim=(/1, 3, (1, dim=3,nvars)/), &
            centering=(/(0, dim=1,nvars)/), &
            varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,nvars) /), &
            vars=(/ C_LOC(grid%primitive(1,1,1,grid%mIdx,0)), C_LOC(velocity), &
                   (C_LOC(grid%primitive(1,1,1,grid%eIdx+dim,0)), dim=0,grid%nVariables-grid%eIdx) /))
      end if
   end subroutine

   subroutine writeVelocityToVTK()
      ! For now we only write rho and v to a VTK file
      ! Otherwise there may be lots and lots of variables
      integer :: mesh_dims(3), dim, iVariance
      character(len=16), dimension(max(6,grid%nVariables)), target :: varnames

      character(len=nMaxCharacters), target :: filename
      
      real(wp), dimension(:,:,:,:), allocatable, target :: velocity
      integer :: nvars

      filename = trim(filenameBase) // ".means.vtk"
      write(*,*) "Writing mean density and velocity to file ", trim(filename)
      filename = trim(filename) // C_NULL_CHAR

      mesh_dims=grid%nCells(1:3)
      if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)

      if(grid%nDimensions==0) then ! Don't plot velocity, plot scalars only

         nvars=grid%nVariables-grid%nDimensions
         varnames(1) = "Density" // C_NULL_CHAR
         varnames(2) = "Temperature" // C_NULL_CHAR
         varnames(3) = "Scalar1" // C_NULL_CHAR
         varnames(4) = "Scalar2" // C_NULL_CHAR
         varnames(5) = "Scalar3" // C_NULL_CHAR
         varnames(6:) = "Scalar" // C_NULL_CHAR
         CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
            ub=0_c_int, dims=mesh_dims+1, &
            x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
            y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
            z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
            nvars=nvars, vardim=(/(1, dim=1,nvars)/), &
            centering=(/(0, dim=1,nvars)/), &
            varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,nvars) /), &
            vars=(/ C_LOC(grid%meanPrimitive(1,1,1,grid%mIdx,0)), &
                   (C_LOC(grid%meanPrimitive(1,1,1,grid%eIdx+dim,0)), dim=0,grid%nVariables-grid%eIdx) /))
      
      
      else ! Save a 3D velocity field for plotting in visit

         allocate(velocity(3, grid%nCells(1), grid%nCells(2), grid%nCells(3)))
         do dim=1, grid%nDimensions
            velocity(dim, :, :, :) = grid%meanPrimitive(:, :, :, grid%jIdx1 + dim - 1 , 0)
         end do   
         velocity(grid%nDimensions+1 : 3, :, :, :) = 0.0_wp

         varnames(1) = "Density" // C_NULL_CHAR
         varnames(2) = "Velocity" // C_NULL_CHAR
         CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
            ub=0_c_int, dims=mesh_dims+1, &
            x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
            y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
            z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
            nvars=2, vardim=(/1, 3/), centering=(/0, 0/), &
            varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1,2) /), &
            vars=(/ C_LOC(grid%meanPrimitive(1,1,1,grid%mIdx,0)), C_LOC(velocity) /))
         
      end if   

   end subroutine

   subroutine writeVariancesToVTK()
      integer :: mesh_dims(3), dim, iVariance
      integer(c_int), dimension(grid%nVariances) :: zeros, ones
      character(len=16), dimension(grid%nVariances), target :: varnames

      character(len=nMaxCharacters), target :: filename

      zeros=0
      ones=1
      filename = trim(filenameBase) // ".variances.vtk"
      write(*,*) "Writing variances to file ", trim(filename)
      filename = trim(filename) // C_NULL_CHAR

      mesh_dims=grid%nCells(1:3)
      if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)

      do iVariance=1, grid%nVariances
         write(varnames(iVariance),"(A,I1,A)") "var", iVariance, C_NULL_CHAR
      end do                  

      CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
         ub=0_c_int, dims=mesh_dims+1, &
         x=(/ (grid%systemLength(1)/grid%nCells(1)*dim, dim=0,mesh_dims(1)) /), &
         y=(/ (grid%systemLength(2)/grid%nCells(2)*dim, dim=0,mesh_dims(2)) /), &
         z=(/ (grid%systemLength(3)/grid%nCells(3)*dim, dim=0,mesh_dims(3)) /), &            
         nvars=grid%nVariances, vardim=ones, centering=zeros, &
         varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1, grid%nVariances) /), &
         vars=(/ (C_LOC(grid%variances(1,1,1,dim)), dim=1, grid%nVariances) /))

   end subroutine

end subroutine  

! This saves the static structure factors to a file
subroutine writeStructureFactors(grid,filenameBase)
   type (HydroGrid), target :: grid
   character(len=*) :: filenameBase

   integer :: structureFactorFile(3) = (/521, 522, 523/)
   integer :: i, j, k, iSpecies, species, species1, species2, nTemps, dim, &
              iVariable, variable, variable1, variable2, iVariance, minky, maxky
   integer, dimension(nMaxDims) ::  mink, maxk, mesh_dims         
   integer :: iCell, jCell, kCell, iDimension, n_fft_cells, &
      iStructureFactor, jStructureFactor, iCross, iFile
   integer, dimension(grid%nStructureFactors) :: varianceIndex
   real(wp) :: wavevector(nMaxDims)

   integer :: ijk(nMaxDims), nModes, iMode, dim1, dim2
   real(wp) :: k2_discrete, modes(nMaxDims,nMaxDims), k2_proj
   real(wp), dimension(nMaxDims) :: kdrh, k_discrete
   
   integer :: n_S_k_vars
   real(wp) :: S_k_av, S_k_variance
   real(wp), allocatable, target :: S_k_array(:,:,:,:)
   character(len=16), dimension(16), target :: varnames
   
   character(25) :: id_string
   character (nMaxCharacters), target :: filename
   
   logical, parameter :: calculateCorrelation=.false.

   if (grid%nStructureFactors<=0) return 
   
   mink=grid%minWavenumber
   maxk=grid%maxWavenumber
   wavevector = 2 * pi / grid%systemLength
   if(.not.grid%periodic) then
      wavevector(2) = grid%systemLength(2) / grid%nCells(2)
      minky=1; maxky=grid%nCells(2)
   else   
      minky=mink(2); maxky=maxk(2)
   end if

   ! Calculate the cross-covariances before we subtract the mean
   EstCovs: if(grid%estimateCovariances) then ! Correlations between different structure factors
      ! This might be useful in variance reduction  
      
      ! Subtract the product of the means to get the covariance:
      iCross=0
      do iStructureFactor = 1, grid%nStructureFactors
      do jStructureFactor = 1, grid%nStructureFactors
      if(iStructureFactor<=jStructureFactor) then
         iCross=iCross+1
         if(iStructureFactor==jStructureFactor) varianceIndex(iStructureFactor) = iCross
         grid%structureFactorsCovariance(:, :, :, iCross) = grid%structureFactorsCovariance(:, :, :, iCross) - &
            real(grid%structureFactors(:, :, :, iStructureFactor)) * real(grid%structureFactors(:, :, :, jStructureFactor))
      end if
      end do
      end do   
      write(*,*) "Variance indexes for structure factors: ", varianceIndex
      
      CalcCorr: if(calculateCorrelation) then ! Calculate correlation coefficients instead of covariances

         ! First the off-diagonal entries   
         iCross=0
         do iStructureFactor = 1, grid%nStructureFactors
         do jStructureFactor = 1, grid%nStructureFactors
         if(iStructureFactor<=jStructureFactor) then
            iCross=iCross+1
            if(iStructureFactor<jStructureFactor) then ! Convert to correlation coefficient         
               grid%structureFactorsCovariance(:, :, :, iCross) = grid%structureFactorsCovariance(:, :, :, iCross) / &
                  sqrt( grid%structureFactorsCovariance(:, :, :, varianceIndex(iStructureFactor)) * &
                        grid%structureFactorsCovariance(:, :, :, varianceIndex(jStructureFactor)) )
            end if   
         end if
         end do
         end do 

         ! Now the diagonal entries      
         if(.false.) then
            do iStructureFactor = 2, grid%nStructureFactors ! Calculate signal/noise ratio |sigma(X)/X|^2
               iCross = varianceIndex(iStructureFactor)
               grid%structureFactorsCovariance(:, :, :, iCross) = &
                  real(grid%structureFactors(:, :, :, iStructureFactor)) * &
                  real(conjg(grid%structureFactors(:, :, :, iStructureFactor))) / &
                  grid%structureFactorsCovariance(:, :, :, iCross)
            end do
         else ! The first pair is special, it is our target variable, all the rest are possible control variables      
            ! Compute the optimal multiplier for all of the control variables:
            do iStructureFactor = 2, grid%nStructureFactors
               iCross = varianceIndex(iStructureFactor)
               grid%structureFactorsCovariance(:, :, :, iCross) =  &
                  grid%structureFactorsCovariance(:, :, :, iStructureFactor) * &
                  sqrt(grid%structureFactorsCovariance(:, :, :, varianceIndex(1)) / &
                       grid%structureFactorsCovariance(:, :, :, iCross))
            end do         
         end if

      else

         ! Normalize by the variance of the special (first) pair:
         do iStructureFactor = 2, grid%nStructureFactors
            grid%structureFactorsCovariance(:,:,:,iStructureFactor) = &
               grid%structureFactorsCovariance(:,:,:,iStructureFactor) / &
               grid%structureFactorsCovariance(:,:,:,1)
         end do      

      end if  CalcCorr      

      ! Signal to noise ratio for the first variable:
      grid%structureFactorsCovariance(:, :, :, 1) = &
         grid%structureFactors(:, :, :, 1) * conjg(grid%structureFactors(:, :, :, 1)) / &
         grid%structureFactorsCovariance(:, :, :, 1)

      ! The k=0 value is problematic so just zero it out:
      ! The k=0 value is problematic:
      iCell = frequencyToArrayIndex(0, grid%nCells(1))
      jCell = frequencyToArrayIndex(0, grid%nCells(2))
      kCell = frequencyToArrayIndex(0, grid%nCells(3))
      grid%structureFactorsCovariance(iCell, jCell, kCell, :) = 0.0_wp

   end if EstCovs
   
   !-----------------------
   ! Now, we subtract the Fourier Transform (FT) of the mean to get a covariance (structure factor)
   ! It makes no sense to do this if there is only one sample as we will get zero
   if( grid%subtractMeanFT .and. (grid%iStep>1)) then
   if(allocated(grid%meanStaticTransforms)) then
      ! We are keep track of the mean of the FT, so we can just subtract it:

      ! Now, we subtract the Fourier Transform (FT) of the mean to get a covariance (structure factor)
      ! S(A,B)=<FT(A)*CONJ(FT(B))>-<FT(A)>*CONJ(<FT(B)>)
      do iStructureFactor = 1, grid%nStructureFactors
         variable1 = grid%structureFactorPairs(1, iStructureFactor)
         variable2 = grid%structureFactorPairs(2, iStructureFactor)

         grid%structureFactors(:, :, :, iStructureFactor) = grid%structureFactors(:, :, :, iStructureFactor) - &
            grid%meanStaticTransforms(:, :, :, variable1) * conjg(grid%meanStaticTransforms(:, :, :, variable2))            
      end do

   else
      ! Take the FT of the mean, instead of the mean of the FT
      ! The way this is written it is destructive since writing to a file changes the stored values
      ! Therefore one must reset all the counters after writing to a file
      
      if(.not.grid%periodic) stop "Non periodic systems not supported if not storing mean of FFTs"

      ! First we need to calculate the Fourier transforms of the means   
      do iVariable = 1, grid%nVariablesToFFT
         species = grid%variablesToFFT(1, iVariable)
         variable = grid%variablesToFFT(2, iVariable)

         grid%staticFFTarray = grid%meanPrimitive(:, :, :, variable, species)
         call FFTW_Execute(grid%staticFFTplan)
         grid%staticFourierTransforms(:, :, :, iVariable) = grid%staticFFTarray

      end do      

      ! S(A,B)=<FT(A)*CONJ(FT(B))>-FT(<A>)*CONJ(FT(<B>))
      do iStructureFactor = 1, grid%nStructureFactors
         variable1 = grid%structureFactorPairs(1, iStructureFactor)
         variable2 = grid%structureFactorPairs(2, iStructureFactor)

         grid%structureFactors(:, :, :, iStructureFactor) = grid%structureFactors(:, :, :, iStructureFactor) - &
            grid%staticFourierTransforms(:, :, :, variable1) * &
            conjg(grid%staticFourierTransforms(:, :, :, variable2))
      end do
   end if
   end if
   
   ! The definition of the theoretical S_k requires normalization of the FFTW result by dV/N_cells
   ! There is also a normalization by n_fft_cells due to FFTWs convention
   if(grid%periodic) then
      n_fft_cells = product(grid%nCells)
   else ! There is no FFT along the y axis
      n_fft_cells = product(grid%nCells) / grid%nCells(2)
   end if
   grid%structureFactors = grid%structureFactors / n_fft_cells * &
      grid%structFactMultiplier * product(grid%systemLength) / product(grid%nCells)
   
   
   !-----------------------
   ! Now actually output the files:

   ForEachStructureFactor: do iStructureFactor = 1, grid%nStructureFactors
     
      write(id_string,"(I25)") iStructureFactor
      
      filename=trim(filenameBase) // ".S_k.pair=" // trim(ADJUSTL(id_string))
      write(*,*) "Writing static structure factor to file ", trim(filename) // ".{Re,Im}.dat"
      open (file = trim(filename) // ".Re.dat", unit=structureFactorFile(1), status = "unknown", action = "write")
      open (file = trim(filename) // ".Im.dat", unit=structureFactorFile(2), status = "unknown", action = "write")

      ! i, j, k simply loop from the negative wavenumbers to the positive wavenumbers.
      ! For example, if nCells(1) = 6, i loops from -2 to 3.
      ! iCell = nCells(1) is the negative frequency corresponding to i = -1.
      ! The next value is the zero component frequency, which is located back at iCell = 1.
      iCell = frequencyToArrayIndex(mink(1), grid%nCells(1))
      jCell = frequencyToArrayIndex(mink(2), grid%nCells(2))
      kCell = frequencyToArrayIndex(mink(3), grid%nCells(3))
      
      WriteS_ky: if ( product(grid%nCells) == grid%nCells(grid%axisToPrint) ) then 
         ! One-dimensional system, easy to do:         
         dim = grid%axisToPrint
         
         do iFile=1,2
            write(structureFactorFile(iFile), "(A,2(A, I0,'-',I0))") "# k, S_k for", &
               " ik=", mink(dim), maxk(dim)
            write(structureFactorFile(iFile), "(A,100(g17.9))") "# dk=",  2 * pi / grid%systemLength(dim)
         end do   
         
         ijk = (/iCell,jCell,kCell/)            
         do j = mink(dim), maxk(dim)
            if((.not.grid%periodic).or.(j/=0)) then
               ! Skip k=0 here since it is easy to do and useful
            
               do iFile=1,2
               if(grid%periodic) then ! Wavenumber
                  write(structureFactorFile(iFile), '(g17.9)', advance="no")  2 * pi * j / grid%systemLength(dim)
               else ! Height
                  write(structureFactorFile(iFile), '(g17.9)', advance="no")  &
                     (j-mink(dim)+0.5_wp) * grid%systemLength(dim) / grid%nCells(dim)
               end if   
               end do

               write(structureFactorFile(1), '(g17.9)', advance="no")  &
                  real(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))
               write(structureFactorFile(2), '(g17.9)', advance="no")  &
                  aimag(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))

               ! New lines:
               write(structureFactorFile(1),*)
               write(structureFactorFile(2),*)
            end if      
            
            ! Find the approriate FFTW frequency index:
            ijk(dim) = ijk(dim) + 1 ; if (ijk(dim) > grid%nCells(dim)) ijk(dim) = 1
            iCell=ijk(1); jCell=ijk(2); kCell=ijk(3); 
         end do

         close (structureFactorFile(1))
         close (structureFactorFile(2))
      
      else if(grid%axisToPrint==2) then ! In higher dimensions we only support special cases
      
         do iFile=1,2
            write(structureFactorFile(iFile), "(A,2(A, I0,'-',I0))") "# ky, S_k for", &
               " ikz=", mink(3), maxk(3), & 
               " ikx=", mink(1), maxk(1)
            write(structureFactorFile(iFile), "(A,100(g17.9))") "# dk=",  2 * pi / grid%systemLength
         end do   
                     
         do j = mink(2), maxk(2)

            do iFile=1,2
               if(grid%periodic) then ! Wavenumber
                  write(structureFactorFile(iFile), '(g17.9)', advance="no")  2 * pi * j / grid%systemLength(2)
               else ! Height
                  write(structureFactorFile(iFile), '(g17.9)', advance="no")  &
                     (j-mink(2)+0.5_wp) * grid%systemLength(2) / grid%nCells(2)
               end if   
            end do   

            do k = mink(3), maxk(3)
            do i = mink(1), maxk(1)
               write(structureFactorFile(1), '(g17.9)', advance="no")  &
                  real(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))
               write(structureFactorFile(2), '(g17.9)', advance="no")  &
                  aimag(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))

            iCell = iCell + 1 ; if (iCell > grid%nCells(1)) iCell = 1
            end do
            kCell = kCell + 1 ; if (kCell > grid%nCells(3)) kCell = 1
            end do

            ! New lines:
            write(structureFactorFile(1),*)
            write(structureFactorFile(2),*)
         
         jCell = jCell + 1 ; if (jCell > grid%nCells(2)) jCell = 1
         end do   

         close (structureFactorFile(1))
         close (structureFactorFile(2))
      
      else WriteS_ky   
      
         !write(0,*) "Only grid%axisToPrint==2 implemented for now when writing S_k.dat files"
         do iFile=1,2
            write(structureFactorFile(iFile), "(A)") &
               "Only grid%axisToPrint==2 implemented for now when writing S_k.dat files"
         end do   

      end if WriteS_ky
   
   end do ForEachStructureFactor   

   if(grid%estimateCovariances) then ! Correlations between different structure factors
      ! This might be useful in variance reduction 
      grid%axisToPrint=2 ! Only y axis implemented here
      
      iCross=0
      do iStructureFactor = 1, grid%nStructureFactors
      do jStructureFactor = 1, grid%nStructureFactors
      if(iStructureFactor<=jStructureFactor) then
         iCross=iCross+1
         
         write(id_string,"(I25)") iCross
         filename=trim(filenameBase) // ".S_k.cov=" // trim(ADJUSTL(id_string))
         write(*,*) "Writing covariances of static structure factors ", iStructureFactor, jStructureFactor , &
            " to file ", trim(filename) // ".dat"
         open (file = trim(filename) // ".dat", unit=structureFactorFile(1), status = "unknown", action = "write")

         do j = mink(2), maxk(2)

         do k = mink(3), maxk(3)
         do i = mink(1), maxk(1)
            write(structureFactorFile(1), '(g17.9)', advance="no")  &
               (grid%structureFactorsCovariance(iCell, jCell, kCell, iCross))

         iCell = iCell + 1 ; if (iCell > grid%nCells(1)) iCell = 1
         end do
         kCell = kCell + 1 ; if (kCell > grid%nCells(3)) kCell = 1
         end do

         ! New lines:
         write(structureFactorFile(1),*)

         jCell = jCell + 1 ; if (jCell > grid%nCells(2)) jCell = 1
         end do   

         close (structureFactorFile(1))  
         
      end if
      end do
      end do      
   end if

   ! Now also write this as VTK files if asked for:
   WriteVTK_S_k: if(grid%writeSpectrumVTK) then   

      ! In order to provide some quantiative values instead of just VTK files for each structure factor we
      ! write the mean and standard deviations over all wavenumbers (excluding the origin)
      filename=trim(filenameBase) // ".S_k.stats.dat"
      write(*,*) "Writing statistics about static structure factors to file ", trim(filename)
      open (file = trim(filename), unit=structureFactorFile(3), status = "unknown", action = "write")
      write(structureFactorFile(3), "(A)") "# Columns are: mean, std, minval, maxval"

      if(writeAbsValue) then ! Only write the absolute value of the static factors  
         n_S_k_vars=1 ! For abs value
         varnames(1)="S_k" // C_NULL_CHAR         
      else
         n_S_k_vars=2 ! One for real and one for imaginary part
         varnames(1)="Real" // C_NULL_CHAR
         varnames(2)="Imag" // C_NULL_CHAR
      end if

      nTemps=n_S_k_vars ! Number of temporary arrays needed
      
      if(writeTheory>=0) then ! Calculate the discrete structure factor theory
         ! Write the actual correct values
         n_S_k_vars = n_S_k_vars + 1
         varnames(n_S_k_vars)="Theory" // C_NULL_CHAR
         nTemps=n_S_k_vars
      else if(any(grid%vectorFactors/=0)) then
         ! We will project the solution onto the discrete modes
         ! Projection onto divergence-free vector fields:
         varnames(n_S_k_vars+1)="Longitudinal" // C_NULL_CHAR
         ! Projection onto the first non-divergence free mode:
         varnames(n_S_k_vars+2)="Vortical-Z" // C_NULL_CHAR
         if(grid%nDimensions>2) then
            ! Projection onto the second non-divergence free mode and the mode-mode correlation:
            varnames(n_S_k_vars+3)="Vortical-XY" // C_NULL_CHAR         
            varnames(n_S_k_vars+4)="Vortical-XY-Z" // C_NULL_CHAR
            nModes=4
         else
            nModes=2      
         end if         
         nTemps=n_S_k_vars+nModes
      end if

      if(any(grid%vectorFactors/=0)) then
         nTemps=nTemps+1 ! Also calculate the trace as an alternative
         varnames(nTemps)="Trace" // C_NULL_CHAR
      end if         
      allocate(S_k_array(mink(1):maxk(1), mink(2):maxk(2), mink(3):maxk(3), nTemps))      

      S_k_array = 0.0_wp ! Initialize the sums      
      do iStructureFactor = 1, grid%nStructureFactors         
         
         ! First we need to copy into a temporary array:
         iCell = frequencyToArrayIndex(mink(1), grid%nCells(1))
         jCell = frequencyToArrayIndex(mink(2), grid%nCells(2))
         kCell = frequencyToArrayIndex(mink(3), grid%nCells(3))
         do k = mink(3), maxk(3)
         do j = mink(2), maxk(2)
         do i = mink(1), maxk(1)
            ijk = (/i,j,k/)
            
            kdrh = ijk * pi / grid%nCells ! This is k*dx/2
            k_discrete = ijk * 2*pi / grid%systemLength ! Wavenumber
            if(abs(writeTheory)==1) then ! Account for discretization artifacts
               where(ijk==0)
                  k_discrete = 0
               else where
                  k_discrete = k_discrete * (sin(kdrh)/kdrh)
               end where
            end if   
            k2_discrete = max(sum(k_discrete**2), epsilon(1.0_wp))
            
            ! If we want to project velocities onto divergence-free modes:         
            k2_proj = max(k2_discrete - k_discrete(3)**2, epsilon(1.0_wp)) ! kx^2+ky^2
            modes(1,:) = k_discrete / sqrt(k2_discrete) ! The divergence of velocity (zero mode)
            if((i==0).and.(j==0)) then ! This is a singular limit since it depends on the angle between k and the x/y axes
               ! First incompressible mode:
               modes(2,:) = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
               ! Second incompressible mode:
               modes(3,:) = (/ 0.0_wp, 1.0_wp, 0.0_wp /)
            else
               ! First incompressible mode:
               modes(2,:) = (/ -k_discrete(2), k_discrete(1), 0.0_wp /) / sqrt(k2_proj)
               ! Second incompressible mode:
               modes(3,:) = (/ -k_discrete(3)*k_discrete(1), -k_discrete(3)*k_discrete(2), k2_proj /) / sqrt(k2_discrete*k2_proj)
            end if
               
            if(writeAbsValue) then
               S_k_array(i,j,k,1)=abs(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))
            else
               S_k_array(i,j,k,1)=real(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))
               S_k_array(i,j,k,2)=aimag(grid%structureFactors(iCell, jCell, kCell, iStructureFactor))
            end if   
            
            if(writeTheory==0) then ! Write the continuum theoretical prediction for incompressible hydro
               if(grid%vectorFactors(iStructureFactor)>0) then ! <v_x, v_x> = P_xx
                  S_k_array(i,j,k,n_S_k_vars) = 1.0_wp-kdrh(iStructureFactor)**2/max(sum(kdrh**2), epsilon(1.0_wp))
               else if(grid%vectorFactors(iStructureFactor)<0) then ! <v_x, v_y> = P_xy
                  S_k_array(i,j,k,n_S_k_vars) = -kdrh(1)*kdrh(2)/max(sum(kdrh**2), epsilon(1.0_wp))
               end if   
            else if(writeTheory>0) then ! Write the discrete theoretical prediction (MAC projection)
               if(grid%vectorFactors(iStructureFactor)>0) then ! <v_x, v_x> = P_xx
                  if((i==0).and.(j==0).and.(k==0)) then
                     S_k_array(i,j,k,n_S_k_vars) = 1.0_wp
                  else   
                     S_k_array(i,j,k,n_S_k_vars) = 1.0_wp-sin(kdrh(iStructureFactor))**2 / max(sum(sin(kdrh)**2), epsilon(1.0_wp))
                  end if
               else if(grid%vectorFactors(iStructureFactor)<0) then ! <v_x, v_y> = P_xy
                  S_k_array(i,j,k,n_S_k_vars) = -sin(kdrh(1))*sin(kdrh(2)) / max(sum(sin(kdrh)**2), epsilon(1.0_wp))
               end if
            else if(any(grid%vectorFactors/=0)) then
               ! Project onto the discretely-divergence free modes
               if(grid%vectorFactors(iStructureFactor)>0) then ! A self-correlation such as <v_x, v_x>
                  dim = grid%vectorFactors(iStructureFactor)
                  ! Add the self contribution from this structure factor to each mode:
                  do iMode=1, grid%nDimensions ! Squared amplitude of a mode
                     S_k_array(i,j,k,n_S_k_vars+iMode) = S_k_array(i,j,k,n_S_k_vars+iMode) + &
                        real((modes(iMode,dim)**2)  * grid%structureFactors(iCell, jCell, kCell, iStructureFactor), wp)
                  end do   
                  if(grid%nDimensions>2) then ! Cross-correlation between two modes (2 and 3)
                     iMode = 4
                     S_k_array(i,j,k,n_S_k_vars+iMode) = S_k_array(i,j,k,n_S_k_vars+iMode) + &
                        real((modes(2,dim)*modes(3,dim)) * &
                           grid%structureFactors(iCell, jCell, kCell, iStructureFactor), wp)                  
                  end if
               else if(grid%vectorFactors(iStructureFactor)<0) then ! A cross-correlation such as <v_x, v_y>
                  select case(grid%vectorFactors(iStructureFactor))
                  case(-1) ! vx-vy
                     dim1=1; dim2=2
                  case(-2) ! vy-vz
                     dim1=2; dim2=3
                  case(-3) ! vx-vz
                     dim1=1; dim2=3
                  case default
                     stop "Negative values of vectorFactors must be between -3 and -1"
                  end select
                  do iMode=1, grid%nDimensions ! Squared amplitude of a mode
                     S_k_array(i,j,k,n_S_k_vars+iMode) = S_k_array(i,j,k,n_S_k_vars+iMode) + &
                        real(2*(modes(iMode,dim1)*modes(iMode,dim2))  * &
                           grid%structureFactors(iCell, jCell, kCell, iStructureFactor), wp)
                  end do
                  if(grid%nDimensions>2) then ! Cross-correlation between two modes (2 and 3)
                     iMode = 4
                     S_k_array(i,j,k,n_S_k_vars+iMode) = S_k_array(i,j,k,n_S_k_vars+iMode) + &
                        real((modes(2,dim1)*modes(3,dim2) + modes(3,dim1)*modes(2,dim2)) * &
                           grid%structureFactors(iCell, jCell, kCell, iStructureFactor), wp)                  
                  end if
               end if
               
            end if            
            
            ! Calculate the trace of the tensor structure factor
            if(grid%vectorFactors(iStructureFactor)>0) then
               S_k_array(i,j,k,nTemps) = S_k_array(i,j,k,nTemps) + S_k_array(i,j,k,1) ! Add the real/abs parts
            end if   

         iCell = iCell + 1 ; if (iCell > grid%nCells(1)) iCell = 1
         end do
         jCell = jCell + 1 ; if (jCell > grid%nCells(2)) jCell = 1
         end do
         kCell = kCell + 1 ; if (kCell > grid%nCells(3)) kCell = 1
         end do
                  
         do dim=1, n_S_k_vars
            call WriteStats_S_k(id=iStructureFactor)
         end do                     

         !-----------------------
         ! Now write the full structure factors to a VTK file:
         
         write(id_string,"(I25)") iStructureFactor
         filename=trim(filenameBase) // ".S_k.pair=" // trim(adjustl(id_string)) // ".vtk"
         write(*,*) "Writing static structure factor to VTK file ", trim(filename)
         filename = trim(filename) // C_NULL_CHAR

         ! We write the real and imaginary part together here:
         mesh_dims=grid%nCells(1:3)
         if(grid%nCells(3)<=1) then
            mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)            
               
            CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
               ub=0_c_int, dims=mesh_dims+1, &
               x=(/ ( wavevector(1)*(dim-0.5_wp), dim = mink(1), maxk(1)+1 ) /), &
               y=(/ ( wavevector(2)*(dim-0.5_wp), dim = minky, maxky+1 ) /), &
               z=(/ 0.0_wp, 0.0_wp /), &
               nvars=n_S_k_vars, vardim=(/(1, dim=1,n_S_k_vars)/), centering=(/(0, dim=1,n_S_k_vars)/), &
               varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1, n_S_k_vars) /), &
               vars=(/ (C_LOC(S_k_array(mink(1),mink(2),mink(3),dim)), dim=1, n_S_k_vars) /))
         
         else if (maxky==minky) then ! Also a 2D result, so we swap z and y in the vtk file for easier rendering
         
            mesh_dims(2)=grid%nCells(3)
            mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK) 

            CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
               ub=0_c_int, dims=mesh_dims+1, &
               x=(/ ( wavevector(1)*(dim-0.5_wp), dim = mink(1), maxk(1)+1 ) /), &
               y=(/ ( wavevector(3)*(dim-0.5_wp), dim = mink(3), maxk(3)+1 ) /), &
               z=(/ 0.0_wp, 0.0_wp /), &
               nvars=n_S_k_vars, vardim=(/(1, dim=1,n_S_k_vars)/), centering=(/(0, dim=1,n_S_k_vars)/), &
               varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1, n_S_k_vars) /), &
               vars=(/ (C_LOC(S_k_array(mink(1),mink(2),mink(3),dim)), dim=1, n_S_k_vars) /))
             
         else

            CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
               ub=0_c_int, dims=mesh_dims+1, &
               x=(/ ( wavevector(1)*(dim-0.5_wp), dim = mink(1), maxk(1)+1 ) /), &
               y=(/ ( wavevector(2)*(dim-0.5_wp), dim = minky, maxky+1 ) /), &
               z=(/ ( wavevector(3)*(dim-0.5_wp), dim = mink(3), maxk(3)+1 ) /), &
               nvars=n_S_k_vars, vardim=(/(1, dim=1,n_S_k_vars)/), centering=(/(0, dim=1,n_S_k_vars)/), &
               varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=1, n_S_k_vars) /), &
               vars=(/ (C_LOC(S_k_array(mink(1),mink(2),mink(3),dim)), dim=1, n_S_k_vars) /))
         
         end if      
               

      end do
      
      if(nTemps > n_S_k_vars) then ! Write the trace separately
      
         write(structureFactorFile(3),*) ! Newline         
         do dim=n_S_k_vars+1,nTemps
            call WriteStats_S_k(id=dim-n_S_k_vars)
         end do   

         filename=trim(filenameBase) // ".S_k.trace.vtk"
         write(*,*) "Writing static structure factor to VTK file ", trim(filename)
         filename = trim(filename) // C_NULL_CHAR

         mesh_dims=grid%nCells(1:3)
         if(grid%nCells(3)<=1) mesh_dims(3)=0 ! Indicate that this is a 2D grid (sorry, no 1D grid in VTK)            
         
         CALL WriteRectilinearVTKMesh(filename=C_LOC(filename(1:1)), &
            ub=0_c_int, dims=mesh_dims+1, &
            x=(/ ( wavevector(1)*(dim-0.5_wp), dim = mink(1), maxk(1)+1 ) /), &
            y=(/ ( wavevector(2)*(dim-0.5_wp), dim = minky, maxky+1 ) /), &
            z=merge(0.0_wp, 1.0_wp, mesh_dims(3)==0) * &
              (/ ( wavevector(3)*(dim-0.5_wp), dim = mink(3), maxk(3)+1 ) /), &
            nvars=nTemps-n_S_k_vars, vardim=(/(1, dim=n_S_k_vars+1,nTemps)/), &
            centering=(/(0, dim=n_S_k_vars+1,nTemps)/), &
            varnames=(/ (C_LOC(varnames(dim)(1:1)), dim=n_S_k_vars+1,nTemps) /), &
            vars=(/ (C_LOC(S_k_array(mink(1),mink(2),mink(3),dim)), dim=n_S_k_vars+1,nTemps) /))
      
      end if      

      close(unit=structureFactorFile(3))
      
   end if WriteVTK_S_k

   if (grid%iTimeSeries > 0) call writeDynamicFactors(grid,filenameBase)

contains

   subroutine WriteStats_S_k(id)
      integer, intent(in) :: id
      
      ! The k=0 wavenumber is not meaningful
      ! Instead, we replace it with the average of the rest of the structure factors

      S_k_av = sum (S_k_array(:, :, :, dim))
      S_k_av = (S_k_av - S_k_array(0,0,0, dim)) / &
         max(1, size(grid%structureFactors(:, :, :, dim)) - 1 )
      S_k_array(0,0,0, dim) = S_k_av

      S_k_variance = sum ((S_k_array(:, :, :, dim)-S_k_av)**2) / &
         max(1, size(grid%structureFactors(:, :, :, dim)) - 1 )

      write(structureFactorFile(3), "(A,I0,2A)") "# id =", id, &
         " var = ", varnames(dim)(1:len_trim(varnames(dim))-1)
      write(structureFactorFile(3), "(100(g17.9))") S_k_av, sqrt(S_k_variance), &
         minval(S_k_array(:, :, :, dim)), maxval(S_k_array(:, :, :, dim)) 
           
   end subroutine  
    
end subroutine writeStructureFactors

! This saves the dynamic structure factors to a file
subroutine writeDynamicFactors(grid,filenameBase)
   type (HydroGrid), target :: grid
   character(len=*) :: filenameBase

   integer :: structureFactorFile(2) = (/ 621, 622 /)
   integer :: i, j, k, iSpecies, species, species1, species2, nTemps, dim, &
              iVariable, variable, variable1, variable2, iVariance
   integer :: iCell, jCell, kCell, wCell, iDimension, iFile, &
      iStructureFactor, jStructureFactor, iWave, iWavenumber
   real(wp) :: wavevector(nMaxDims), wavefrequency
   
   character(25) :: id_string
   character (nMaxCharacters), target :: filename
   
   ! Dynamic factors will also be saved for the vortical (solenoidal) modes and divergence of velocity:
   integer :: ijk(nMaxDims), nModes, iMode, dim1, dim2
   real(wp) :: k2_discrete, modes(nMaxDims,nMaxDims), k2_proj
   real(wp), dimension(nMaxDims) :: kdrh, k_discrete
   complex(wp) :: velocityModes(nMaxDims, grid%nSavedSnapshots)
   
   if (grid%nWavenumbers <= 0) return
         
   ! Recall grid%dynamicFactors(grid%nSavedSnapshots, grid%nWavenumbers, grid%nStructureFactors)

   ! The definition of the theoretical S_k_w requires normalization of the FFTW result by dx/N_cells:
   grid%dynamicFactors = grid%dynamicFactors * grid%structFactMultiplier * &
      product(grid%systemLength) / real(product(grid%nCells),wp)**2 * &
      grid%timestep / grid%nSavedSnapshots ! Note that product(grid%nCells)^2 can overflow in integer
      
   ! For projecting the velocity onto solenoidal modes:
   nModes=grid%nDimensions ! We do not compute cross-correlations here
   
   do iWavenumber = 1, grid%nWavenumbers
   
      write(id_string,"(I25)") iWavenumber
      filename=trim(filenameBase) // ".S_k_w.k=" // trim(ADJUSTL(id_string))
      
      write(*,*) "Writing dynamic structure factor S(k,w) for k=", grid%selectedWavenumbers(:, iWavenumber), &
         " to file ", trim(filename) // ".{Re,Im}.dat"
      open (file = trim(filename) // ".Re.dat", unit=structureFactorFile(1), status = "unknown", action = "write")
      open (file = trim(filename) // ".Im.dat", unit=structureFactorFile(2), status = "unknown", action = "write")
      
      ! This only works if the selected wavenumbers are positive, but this is always the case anyway:
      wavevector = 2 * pi * (grid%selectedWavenumbers(:, iWavenumber)) / grid%systemLength
      if(.not.grid%periodic) then ! This is actually y coordinate not k_y
         wavevector(2) = (grid%selectedWavenumbers(2, iWavenumber)+0.5_wp) * grid%systemLength(dim) / grid%nCells(dim)
      end if
      do iFile=1,2
         write(structureFactorFile(iFile), '(A,100G17.9)')  "# k=", wavevector            
      end do   

      ProjectModes: if(any(grid%vectorFactors/=0)) then
         ! We also compute vortical and solenoidal modes for velocity
         
         ijk = grid%selectedWavenumbers(:, iWavenumber) ! Wave-index
         
         kdrh = ijk * pi / grid%nCells ! This is k*dx/2
         k_discrete = ijk * 2*pi / grid%systemLength ! Wavenumber
         if(abs(writeTheory)==1) then ! Account for discretization artifacts
            where(ijk==0)
               k_discrete = 0
            else where
               ! Consider discrete divergence-free vector fields:
               k_discrete = k_discrete * (sin(kdrh)/kdrh)
            end where
         end if
         k2_discrete = max(sum(k_discrete**2), epsilon(1.0_wp))

         ! If we want to project velocities onto divergence-free modes:         
         k2_proj = max(k2_discrete - k_discrete(3)**2, epsilon(1.0_wp)) ! kx^2+ky^2
         modes(1,:) = k_discrete / sqrt(k2_discrete) ! The divergence of velocity (zero mode)
         if(all(ijk(1:2)==0)) then ! This is a singular limit since it depends on the angle between k and the x/y axes
            ! First incompressible mode:
            modes(2,:) = (/ 1.0_wp, 0.0_wp, 0.0_wp /)
            ! Second incompressible mode:
            modes(3,:) = (/ 0.0_wp, 1.0_wp, 0.0_wp /)
         else
            ! First incompressible mode:
            modes(2,:) = (/ -k_discrete(2), k_discrete(1), 0.0_wp /) / sqrt(k2_proj)
            ! Second incompressible mode:
            modes(3,:) = (/ -k_discrete(3)*k_discrete(1), -k_discrete(3)*k_discrete(2), k2_proj /) / sqrt(k2_discrete*k2_proj)
         end if
         !write(*,*) "ijk=", ijk, " div=", modes(1,1:grid%nDimensions), " vort=", modes(2,1:grid%nDimensions)
         
         velocityModes = 0.0_wp
         do iStructureFactor = 1, grid%nStructureFactors               
            ! Project onto the discretely-divergence free modes into velocityModes(nMaxDims,grid%nSavedSnapshots)
            
            ! This assumes the dynamics is time reversible so omega and -omega are equivalent            
            if(grid%vectorFactors(iStructureFactor)>0) then ! A self-correlation such as <v_x, v_x>
               dim = grid%vectorFactors(iStructureFactor)
               ! Add the self contribution from this structure factor to each mode:
               do iMode=1, grid%nDimensions ! Squared amplitude of a mode
                  velocityModes(iMode,:) = velocityModes(iMode,:) + &
                     (modes(iMode,dim)**2)  *  grid%dynamicFactors(:, iWavenumber, iStructureFactor) 
                  !if(iMode==2) write(*,*) iMode, " Adding transform ", iStructureFactor, &
                  !   " with coefficient ", (modes(iMode,dim)**2)
               end do   
            else if(grid%vectorFactors(iStructureFactor)<0) then ! A cross-correlation such as <v_x, v_y>
               select case(grid%vectorFactors(iStructureFactor))
               case(-1) ! vx-vy
                  dim1=1; dim2=2
               case(-2) ! vy-vz
                  dim1=2; dim2=3
               case(-3) ! vx-vz
                  dim1=1; dim2=3
               case default
                  stop "Negative values of vectorFactors must be between -3 and -1"
               end select
               do iMode=1, grid%nDimensions ! Squared amplitude of a mode
                  velocityModes(iMode,:) = velocityModes(iMode,:) + &
                     modes(iMode,dim1)*modes(iMode,dim2) * grid%dynamicFactors(:, iWavenumber, iStructureFactor) + &
                     modes(iMode,dim1)*modes(iMode,dim2) * conjg(grid%dynamicFactors(:, iWavenumber, iStructureFactor))
                  !if(iMode==2) write(*,*) iMode, " Adding transform ", iStructureFactor, &
                  !   " with coefficient ", 2*modes(iMode,dim1)*modes(iMode,dim2) 
               end do
            end if
            
         end do
         
      end if ProjectModes

      wCell = frequencyToArrayIndex(grid%minWavefrequency, grid%nSavedSnapshots)
      do iWave=grid%minWavefrequency, grid%maxWavefrequency         
         wavefrequency = 2 * pi * iWave / (grid%nSavedSnapshots * grid%timestep)
         
         if(iWave/=0) then ! There is no point in writing zero frequency here
            write(structureFactorFile(1), '(1000g17.9)', advance="no") wavefrequency, &
               real(grid%dynamicFactors(wCell, iWavenumber, :))
            write(structureFactorFile(2), '(1000g17.9)', advance="no") wavefrequency, &
               aimag(grid%dynamicFactors(wCell, iWavenumber, :))
            if(any(grid%vectorFactors/=0)) then
               write(structureFactorFile(1), '(1000g17.9)') real(velocityModes(1:nModes,wCell))
               write(structureFactorFile(2), '(1000g17.9)') aimag(velocityModes(1:nModes,wCell))
            else ! Write end of lines
               write(structureFactorFile(1), '(1000g17.9)')
               write(structureFactorFile(2), '(1000g17.9)')                
            end if   
         end if      
      
         wCell = wCell + 1 ; if (wCell > grid%nSavedSnapshots) wCell = 1

      end do      
        
      close (structureFactorFile(1))
      close (structureFactorFile(2))
   end do    
   
   ! Now also write the correlation functions in real time
   !-----------------------------
   
   ! Unfo the normalization:
   grid%dynamicFactors = grid%dynamicFactors / grid%timestep / grid%nSavedSnapshots
   
   ! Temporarily set the plan to inverse transform:
   call FFTW_PlanDFT(grid%dynamicFFTplan, grid%nSavedSnapshots, &
      grid%dynamicFFTarray, grid%dynamicFFTarray, FFTW_BACKWARD, FFTW_ESTIMATE)
   
   do iWavenumber = 1, grid%nWavenumbers

      do iStructureFactor = 1, grid%nStructureFactors
         grid%dynamicFFTarray = grid%dynamicFactors(:, iWavenumber, iStructureFactor)
         call FFTW_Execute(grid%dynamicFFTplan)
         ! These will be reset anyway, so we can use them as temporary storage:
         grid%dynamicFactors(:, iWavenumber, iStructureFactor) = grid%dynamicFFTarray
      end do
   
      write(id_string,"(I25)") iWavenumber
      filename=trim(filenameBase) // ".S_k_t.k=" // trim(ADJUSTL(id_string))
      
      write(*,*) "Writing dynamic structure factor S(k,t) for k=", grid%selectedWavenumbers(:, iWavenumber), &
         " to file ", trim(filename) // ".{Re,Im}.dat"
      open (file = trim(filename) // ".Re.dat", unit=structureFactorFile(1), status = "unknown", action = "write")
      open (file = trim(filename) // ".Im.dat", unit=structureFactorFile(2), status = "unknown", action = "write")
      
      ! This only works if the selected wavenumbers are positive, but this is always the case anyway:
      do iFile=1,2
         write(structureFactorFile(iFile), '(A,100G17.9)')  "# k=", &
            2 * pi * (grid%selectedWavenumbers(:, iWavenumber)) / grid%systemLength
      end do   

      wCell = frequencyToArrayIndex(grid%minWavefrequency, grid%nSavedSnapshots)
      do iWave=grid%minWavefrequency, grid%maxWavefrequency         
         
         write(structureFactorFile(1), '(1000g17.9)') iWave * grid%timestep, &
            real(grid%dynamicFactors(wCell, iWavenumber, :))
         write(structureFactorFile(2), '(1000g17.9)') iWave * grid%timestep, &
            aimag(grid%dynamicFactors(wCell, iWavenumber, :))
      
         wCell = wCell + 1 ; if (wCell > grid%nSavedSnapshots) wCell = 1
      end do      
      
      close (structureFactorFile(1))
      close (structureFactorFile(2))
   end do

   ! Reset the plan back to forward transform
   call FFTW_PlanDFT(grid%dynamicFFTplan, grid%nSavedSnapshots, &
      grid%dynamicFFTarray, grid%dynamicFFTarray, FFTW_FORWARD, FFTW_ESTIMATE)
      
end subroutine

!------------------------------------------------------------------------
! Utilities:
!------------------------------------------------------------------------      

subroutine energyToTemperature(grid, density, velocity, energy, temperature)
   type (HydroGrid), target :: grid
   real(wp), intent(in), dimension(0:grid%nSpecies-1) :: density
   real(wp), intent(in), dimension(grid%nDimensions, 0:grid%nFluids) :: velocity
   real(wp), intent(in), dimension(0:grid%nFluids) :: energy
   real(wp), intent(out), dimension(0:grid%nFluids)  :: temperature

   real(wp) :: densityOfSpeciesN, totalHeatCapacity
   real(wp), dimension(0:grid%nFluids) :: internalEnergy

   densityOfSpeciesN = density(0) - sum(density(1:))
   totalHeatCapacity = ( sum(grid%heatCapacity(1 :grid%nSpecies-1) * density(1:) ) + &
      grid%heatCapacity(grid%nSpecies) * densityOfSpeciesN ) / density(0) 
   internalEnergy = energy - 0.5_wp * density(0:grid%nFluids) * sum( velocity**2, 1) 

   temperature(0) = internalEnergy(0) / (totalHeatCapacity * density(0) )
   where(density(1:grid%nFluids)>0)
      temperature(1:) = internalEnergy(1:grid%nFluids) / (grid%heatCapacity(1:grid%nFluids) * density(1:grid%nFluids) )
   else where
      temperature(1:) = temperature(0) ! Pretend it is the same as the average temperature
   end where   
   
   !write(*,*) "e2T:", density(0), totalHeatCapacity, temperature(0)
end subroutine energyToTemperature



integer function frequencyToArrayIndex(frequency, nCells) result(arrayIndex)
   integer, intent(in) :: frequency, nCells

   if (frequency >= 0) then
      arrayIndex = frequency + 1
   else
      arrayIndex = nCells + frequency + 1
   end if

   if (arrayIndex > nCells) arrayIndex = 1

end function frequencyToArrayIndex

      
end module HydroGridModule
