MODULE HydroGrid_Interface ! Interface to my HydroGrid module
   ! This one is a *fake* interface that does nothing, so as to avoid having to link in libraries and overheads etc.
   IMPLICIT NONE
   PUBLIC
   
   INTEGER, PARAMETER :: dp = KIND(0.0d0)

CONTAINS

   SUBROUTINE InitializeHydroGrid (n_dims, fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, &
      equilibrium_state, hydro_c_v, render)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT(IN) :: fileunit
      REAL (dp), INTENT (IN) :: box_lengths (n_dims)
      INTEGER, DIMENSION (n_dims), INTENT (IN) :: n_cells
      REAL (dp), INTENT (IN) :: extra_scale
      REAL (dp), INTENT (IN) :: hydro_dt
      LOGICAL, DIMENSION (n_dims), INTENT (IN) :: periodic
      REAL (dp), INTENT (IN) :: wall_speed, plug_speed
      REAL (dp), DIMENSION (2+n_dims), INTENT (IN) :: equilibrium_state
      REAL (dp), INTENT (IN) :: hydro_c_v
      LOGICAL, INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE InitializeHydroGridMixture (n_dims, fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, wall_concentrations, &
      equilibrium_state, hydro_c_v, render)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT(IN) :: fileunit
      REAL (dp), INTENT (IN) :: box_lengths (n_dims)
      INTEGER, DIMENSION (n_dims), INTENT (IN) :: n_cells
      REAL (dp), INTENT (IN) :: extra_scale
      REAL (dp), INTENT (IN) :: hydro_dt
      LOGICAL, DIMENSION (n_dims), INTENT (IN) :: periodic
      REAL (dp), INTENT (IN) :: wall_speed, plug_speed, wall_concentrations (2)
      REAL (dp), DIMENSION (3+n_dims), INTENT (IN) :: equilibrium_state
      REAL (dp), INTENT (IN) :: hydro_c_v (2)
      LOGICAL, INTENT(IN) :: render      
   END SUBROUTINE

   SUBROUTINE UpdateHydroGrid(n_dims, n_hydro_cells, n_momentum_dims, densities, momenta_densities, energy_densities, &
      second_densities, second_momenta)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT (IN) :: n_hydro_cells, n_momentum_dims
      REAL (dp), INTENT (IN) :: densities (n_hydro_cells), momenta_densities (n_momentum_dims,n_hydro_cells), &
         energy_densities (n_hydro_cells)
      REAL (dp), INTENT (IN), OPTIONAL :: second_densities (n_hydro_cells), &
         second_momenta (n_momentum_dims,n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE DestroyHydroGrid (n_dims)
      INTEGER, INTENT(IN) :: n_dims
   END SUBROUTINE

   SUBROUTINE ResetHydroGrid (n_dims)
      INTEGER, INTENT(IN) :: n_dims
   END SUBROUTINE

   SUBROUTINE SaveHydroGrid (n_dims,id)
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT(IN), OPTIONAL :: id
   END SUBROUTINE

   SUBROUTINE TestHydroGrid (n_dims, active)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, INTENT(IN) :: n_dims
      LOGICAL, INTENT(OUT) :: active      
      active=.TRUE.
   END SUBROUTINE

END MODULE
