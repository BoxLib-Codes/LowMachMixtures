MODULE HydroGrid_Interface ! Interface to my HydroGrid module
   USE, INTRINSIC :: ISO_C_BINDING
   !USE HydroGrid_1
   !USE HydroGrid_2
   !USE HydroGrid_3
   IMPLICIT NONE
   PUBLIC

   INTEGER, PARAMETER :: dp = KIND(0.0d0)

INTERFACE   

   SUBROUTINE InitializeHydroGridMixture_1 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, wall_concentrations, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGridMixture_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (1)
      INTEGER (KIND=C_INT), DIMENSION (1), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (1), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed, wall_concentrations (2)
      REAL (KIND=C_DOUBLE), DIMENSION (1+3), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v (2)
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE InitializeHydroGrid_1 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGrid_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (1)
      INTEGER (KIND=C_INT), DIMENSION (1), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (1), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed
      REAL (KIND=C_DOUBLE), DIMENSION (1+2), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridMixture_1 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities, second_momenta) &
      BIND (C, NAME='UpdateHydroGridMixture_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), second_momenta (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridBinary_1 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities) &
      BIND (C, NAME='UpdateHydroGridBinary_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGrid_1 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities) &
      BIND (C, NAME='UpdateHydroGrid_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), momenta_densities (n_dims,n_hydro_cells), &
         energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE DestroyHydroGrid_1 () BIND (C, NAME='DestroyHydroGrid_1')
   END SUBROUTINE
   SUBROUTINE SaveHydroGrid_1 (id) BIND (C, NAME='SaveHydroGrid_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: id
   END SUBROUTINE
   SUBROUTINE ResetHydroGrid_1 () BIND (C, NAME='ResetHydroGrid_1')
   END SUBROUTINE
   SUBROUTINE TestHydroGrid_1 (active) BIND (C, NAME='TestHydroGrid_1')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(KIND=C_INT), INTENT(OUT) :: active
   END SUBROUTINE

!--------

   SUBROUTINE InitializeHydroGridMixture_2 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, wall_concentrations, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGridMixture_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (2)
      INTEGER (KIND=C_INT), DIMENSION (2), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (2), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed, wall_concentrations (2)
      REAL (KIND=C_DOUBLE), DIMENSION (2+3), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v (2)
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE InitializeHydroGrid_2 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGrid_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (2)
      INTEGER (KIND=C_INT), DIMENSION (2), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (2), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed
      REAL (KIND=C_DOUBLE), DIMENSION (2+2), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridMixture_2 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities, second_momenta) &
      BIND (C, NAME='UpdateHydroGridMixture_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), second_momenta (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridBinary_2 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities) &
      BIND (C, NAME='UpdateHydroGridBinary_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGrid_2 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities) &
      BIND (C, NAME='UpdateHydroGrid_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), momenta_densities (n_dims,n_hydro_cells), &
         energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE DestroyHydroGrid_2 () BIND (C, NAME='DestroyHydroGrid_2')
   END SUBROUTINE
   SUBROUTINE SaveHydroGrid_2 (id) BIND (C, NAME='SaveHydroGrid_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: id
   END SUBROUTINE
   SUBROUTINE ResetHydroGrid_2 () BIND (C, NAME='ResetHydroGrid_2')
   END SUBROUTINE
   SUBROUTINE TestHydroGrid_2 (active) BIND (C, NAME='TestHydroGrid_2')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(KIND=C_INT), INTENT(OUT) :: active
   END SUBROUTINE

!--------

   SUBROUTINE InitializeHydroGridMixture_3 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, wall_concentrations, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGridMixture_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (2)
      INTEGER (KIND=C_INT), DIMENSION (3), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (3), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed, wall_concentrations (2)
      REAL (KIND=C_DOUBLE), DIMENSION (3+3), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v (2)
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE InitializeHydroGrid_3 (fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
      periodic, wall_speed, plug_speed, &
      equilibrium_state, hydro_c_v, render) BIND (C, NAME='InitializeHydroGrid_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (C_INT), INTENT(IN) :: fileunit
      REAL (KIND=C_DOUBLE), INTENT (IN) :: box_lengths (2)
      INTEGER (KIND=C_INT), DIMENSION (3), INTENT (IN) :: n_cells
      REAL (KIND=C_DOUBLE), INTENT (IN) :: extra_scale
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_dt
      INTEGER(C_INT), DIMENSION (3), INTENT (IN) :: periodic
      REAL (KIND=C_DOUBLE), INTENT (IN) :: wall_speed, plug_speed
      REAL (KIND=C_DOUBLE), DIMENSION (3+2), INTENT (IN) :: equilibrium_state
      REAL (KIND=C_DOUBLE), INTENT (IN) :: hydro_c_v
      INTEGER(C_INT), INTENT(IN) :: render
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridMixture_3 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities, second_momenta) &
      BIND (C, NAME='UpdateHydroGridMixture_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), second_momenta (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGridBinary_3 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities, &
      second_densities) &
      BIND (C, NAME='UpdateHydroGridBinary_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), second_densities (n_hydro_cells), &
         momenta_densities (n_dims,n_hydro_cells), energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE UpdateHydroGrid_3 (n_hydro_cells, n_dims, densities, momenta_densities, energy_densities) &
      BIND (C, NAME='UpdateHydroGrid_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER (KIND=C_INT), INTENT (IN) :: n_hydro_cells, n_dims
      REAL (KIND=C_DOUBLE), INTENT (IN) :: densities (n_hydro_cells), momenta_densities (n_dims,n_hydro_cells), &
         energy_densities (n_hydro_cells)
   END SUBROUTINE

   SUBROUTINE DestroyHydroGrid_3 () BIND (C, NAME='DestroyHydroGrid_3')
   END SUBROUTINE
   SUBROUTINE SaveHydroGrid_3 (id) BIND (C, NAME='SaveHydroGrid_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(C_INT), INTENT(IN) :: id
   END SUBROUTINE
   SUBROUTINE ResetHydroGrid_3 () BIND (C, NAME='ResetHydroGrid_3')
   END SUBROUTINE
   SUBROUTINE TestHydroGrid_3 (active) BIND (C, NAME='TestHydroGrid_3')
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER(KIND=C_INT), INTENT(OUT) :: active
   END SUBROUTINE

END INTERFACE

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
      
      INTEGER(C_INT) :: periodic_(n_dims), render_
      
      ! Convert the logicals to integers:
      periodic_=MERGE(1,0,periodic)
      render_=MERGE(1,0,render)
      
      SELECT CASE(n_dims)
      CASE(1)
         CALL InitializeHydroGrid_1(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, equilibrium_state, hydro_c_v, render_)
      CASE(2)
         CALL InitializeHydroGrid_2(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, equilibrium_state, hydro_c_v, render_)      
      CASE(3)
         CALL InitializeHydroGrid_3(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, equilibrium_state, hydro_c_v, render_)
      END SELECT      
      
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
      
      INTEGER(C_INT) :: periodic_(n_dims), render_
      
      ! Convert the logicals to integers:
      periodic_=MERGE(1,0,periodic)
      render_=MERGE(1,0,render)
      
      SELECT CASE(n_dims)
      CASE(1)
         CALL InitializeHydroGridMixture_1(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, wall_concentrations, equilibrium_state, hydro_c_v, render_)
      CASE(2)
         CALL InitializeHydroGridMixture_2(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, wall_concentrations, equilibrium_state, hydro_c_v, render_)      
      CASE(3)
         CALL InitializeHydroGridMixture_3(fileunit, box_lengths, n_cells, extra_scale, hydro_dt, &
                 periodic_, wall_speed, plug_speed, wall_concentrations, equilibrium_state, hydro_c_v, render_)
      END SELECT      
      
   END SUBROUTINE

   SUBROUTINE UpdateHydroGrid(n_dims, n_hydro_cells, n_momentum_dims, densities, momenta_densities, energy_densities, &
      second_densities, second_momenta)
      USE, INTRINSIC :: ISO_C_BINDING
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT (IN) :: n_hydro_cells, n_momentum_dims
      REAL (dp), INTENT (IN) :: densities (n_hydro_cells), momenta_densities (n_momentum_dims,n_hydro_cells), &
         energy_densities (n_hydro_cells)
      REAL (dp), INTENT (IN), OPTIONAL :: second_densities (n_hydro_cells), second_momenta (n_momentum_dims,n_hydro_cells)

      IF(PRESENT(second_momenta).AND.PRESENT(second_densities)) THEN
         SELECT CASE(n_dims)
         CASE(1)
            CALL UpdateHydroGridMixture_1(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities, second_momenta)
         CASE(2)
            CALL UpdateHydroGridMixture_2(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities, second_momenta)
         CASE(3)
            CALL UpdateHydroGridMixture_3(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities, second_momenta)
         END SELECT
      ELSE IF(PRESENT(second_densities)) THEN
         SELECT CASE(n_dims)
         CASE(1)
            CALL UpdateHydroGridBinary_1(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities)
         CASE(2)
            CALL UpdateHydroGridBinary_2(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities)
         CASE(3)
            CALL UpdateHydroGridBinary_3(n_hydro_cells, n_momentum_dims, densities, momenta_densities, &
               energy_densities, second_densities)
         END SELECT
      ELSE
         SELECT CASE(n_dims)
         CASE(1)
            CALL UpdateHydroGrid_1(n_hydro_cells, n_momentum_dims, densities, momenta_densities, energy_densities)
         CASE(2)
            CALL UpdateHydroGrid_2(n_hydro_cells, n_momentum_dims, densities, momenta_densities, energy_densities)
         CASE(3)
            CALL UpdateHydroGrid_3(n_hydro_cells, n_momentum_dims, densities, momenta_densities, energy_densities)
         END SELECT      
      END IF               

   END SUBROUTINE

   SUBROUTINE DestroyHydroGrid (n_dims)
      INTEGER, INTENT(IN) :: n_dims
      SELECT CASE(n_dims)
      CASE(1)
         CALL DestroyHydroGrid_1()
      CASE(2)
         CALL DestroyHydroGrid_2()
      CASE(3)
         CALL DestroyHydroGrid_3()         
      END SELECT
   END SUBROUTINE

   SUBROUTINE ResetHydroGrid (n_dims)
      INTEGER, INTENT(IN) :: n_dims
      SELECT CASE(n_dims)
      CASE(1)
         CALL ResetHydroGrid_1()
      CASE(2)
         CALL ResetHydroGrid_2()
      CASE(3)
         CALL ResetHydroGrid_3()         
      END SELECT
   END SUBROUTINE

   SUBROUTINE SaveHydroGrid (n_dims, id)
      INTEGER, INTENT(IN) :: n_dims
      INTEGER, INTENT(IN), OPTIONAL :: id
      
      INTEGER(C_INT) :: id_
      IF(PRESENT(id)) THEN
         id_=id
      ELSE
         id_=-1
      END IF
      
      SELECT CASE(n_dims)
      CASE(1)
         CALL SaveHydroGrid_1(id_)
      CASE(2)
         CALL SaveHydroGrid_2(id_)
      CASE(3)
         CALL SaveHydroGrid_3(id_)    
      END SELECT
   END SUBROUTINE

   SUBROUTINE TestHydroGrid (n_dims, active)
      INTEGER, INTENT(IN) :: n_dims
      LOGICAL, INTENT(OUT) :: active
      
      INTEGER :: active_
      
      SELECT CASE(n_dims)
      CASE(1)
         CALL TestHydroGrid_1(active_)
      CASE(2)
         CALL TestHydroGrid_2(active_)
      CASE(3)
         CALL TestHydroGrid_3(active_)
      END SELECT
      active=(active_>0) ! Convert to Fortran logical
      
   END SUBROUTINE

END MODULE
