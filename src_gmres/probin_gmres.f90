! This module stores the runtime parameters.  The probin_gmres_init() routine is
! used to initialize the runtime parameters

module probin_gmres_module

  use bl_types
  use bl_space

  implicit none

  ! For comments and instructions on how to set the input parameters see namelist section below
  !------------------------------------------------------------- 
  integer   , save :: precon_type,visc_schur_approx
  integer   , save :: mg_verbose,cg_verbose,mg_max_vcycles,mg_minwidth
  integer   , save :: mg_bottom_solver,mg_nsmooths_down,mg_nsmooths_up,mg_nsmooths_bottom
  integer   , save :: mg_max_bottom_nlevels,stag_mg_verbosity,stag_mg_max_vcycles
  integer   , save :: stag_mg_minwidth,stag_mg_bottom_solver,stag_mg_nsmooths_down
  integer   , save :: stag_mg_nsmooths_up,stag_mg_nsmooths_bottom,stag_mg_smoother
  integer   , save :: stag_mg_max_bottom_nlevels,gmres_verbose,gmres_max_outer,gmres_max_inner
  integer   , save :: gmres_max_iter,gmres_min_iter
  real(dp_t), save :: p_norm_weight,scale_factor,mg_rel_tol,stag_mg_omega,stag_mg_rel_tol
  real(dp_t), save :: gmres_rel_tol,gmres_abs_tol

  !------------------------------------------------------------- 
  ! Input parameters controlled via namelist input, with comments
  !------------------------------------------------------------- 

  ! preconditioner type
  ! 1 = projection preconditioner
  !-1 = projection preconditioner with expensive pressure update
  ! 2 = lower triangular preconditioner
  !-2 = lower triangular preconditioner with negative sign
  ! 3 = upper triangular preconditioner
  !-3 = upper triangular preconditioner with negative sign
  ! 4 = Block diagonal preconditioner
  !-4 = Block diagonal preconditioner with negative sign
  namelist /probin_gmres/ precon_type

  ! use the viscosity-based BFBt Schur complement (from Georg Stadler)
  namelist /probin_gmres/ visc_schur_approx

  ! weighting of pressure when computing norms and inner products
  namelist /probin_gmres/ p_norm_weight

  ! scale theta_alpha, beta, gamma, and b_u by this, and then scale x_p by the inverse
  namelist /probin_gmres/ scale_factor

  ! MAC projection solver parameters:
  namelist /probin_gmres/ mg_verbose            ! multigrid verbosity
  namelist /probin_gmres/ cg_verbose            ! BiCGStab (mg_bottom_solver=1) verbosity
  namelist /probin_gmres/ mg_max_vcycles        ! maximum number of V-cycles
  namelist /probin_gmres/ mg_minwidth           ! length of box at coarsest multigrid level
  namelist /probin_gmres/ mg_bottom_solver      ! bottom solver type
                                                ! 0 = smooths only, controlled by mg_nsmooths_bottom
                                                ! 1 = BiCGStab
                                                ! 4 = Fancy bottom solve that coarsens down additionally
                                                !     and then applies mg_nsmooths_bottom smooths
  namelist /probin_gmres/ mg_nsmooths_down      ! number of smooths at each level on the way down
  namelist /probin_gmres/ mg_nsmooths_up        ! number of smooths at each level on the way up
  namelist /probin_gmres/ mg_nsmooths_bottom    ! number of smooths at the bottom (only if mg_bottom_solver=0)
  namelist /probin_gmres/ mg_max_bottom_nlevels ! for mg_bottom_solver=4, number of additional levels of multigrid
  namelist /probin_gmres/ mg_rel_tol            ! relative tolerance stopping criteria

  ! Staggered multigrid solver parameters
  namelist /probin_gmres/ stag_mg_verbosity       ! verbosity
  namelist /probin_gmres/ stag_mg_max_vcycles     ! max number of v-cycles
  namelist /probin_gmres/ stag_mg_minwidth        ! length of box at coarsest multigrid level
  namelist /probin_gmres/ stag_mg_bottom_solver   ! bottom solver type
                                                  ! 0 = smooths only, controlled by mg_nsmooths_bottom
                                                  ! 4 = Fancy bottom solve that coarsens additionally
                                                  !     and then applies stag_mg_nsmooths_bottom smooths
  namelist /probin_gmres/ stag_mg_nsmooths_down   ! number of smooths at each level on the way down
  namelist /probin_gmres/ stag_mg_nsmooths_up     ! number of smooths at each level on the way up
  namelist /probin_gmres/ stag_mg_nsmooths_bottom ! number of smooths at the bottom
  namelist /probin_gmres/ stag_mg_max_bottom_nlevels ! for stag_mg_bottom_solver=4, number of additional levels of multigrid
  namelist /probin_gmres/ stag_mg_omega           ! weighted-jacobi omega coefficient
  namelist /probin_gmres/ stag_mg_smoother        ! 0 = jacobi; 1 = 2*dm-color Gauss-Seidel
  namelist /probin_gmres/ stag_mg_rel_tol         ! relative tolerance stopping criteria

  ! GMRES solver parameters
  namelist /probin_gmres/ gmres_rel_tol         ! relative tolerance stopping criteria
  namelist /probin_gmres/ gmres_abs_tol         ! absolute tolerance stopping criteria
  namelist /probin_gmres/ gmres_verbose         ! gmres verbosity; if greater than 1, more residuals will be printed out
  namelist /probin_gmres/ gmres_max_outer       ! max number of outer iterations
  namelist /probin_gmres/ gmres_max_inner       ! max number of inner iterations, or restart number
  namelist /probin_gmres/ gmres_max_iter        ! max number of gmres iterations
  namelist /probin_gmres/ gmres_min_iter        ! min number of gmres iterations
  !------------------------------------------------------------- 

contains

  subroutine probin_gmres_init()

    use f2kcli
    use parallel
    use bc_module
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer    :: narg, farg

    character(len=128) :: fname
    character(len=128) :: probin_env

    logical :: lexist
    logical :: need_inputs

    integer :: un, ierr

    narg = command_argument_count()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! initialize the runtime parameters
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Defaults
    precon_type = 1

    visc_schur_approx = 0

    p_norm_weight = 1.d0

    scale_factor = 1.d0

    mg_verbose = 0
    cg_verbose = 0
    mg_max_vcycles = 1
    mg_minwidth = 2
    mg_bottom_solver = 4
    mg_nsmooths_down = 2
    mg_nsmooths_up = 2
    mg_nsmooths_bottom = 8
    mg_max_bottom_nlevels = 10
    mg_rel_tol = 1.d-9

    stag_mg_verbosity = 0
    stag_mg_max_vcycles = 1
    stag_mg_minwidth = 2
    stag_mg_bottom_solver = 4
    stag_mg_nsmooths_down = 2
    stag_mg_nsmooths_up = 2
    stag_mg_nsmooths_bottom = 8
    stag_mg_max_bottom_nlevels = 10
    stag_mg_omega = 1.d0
    stag_mg_smoother = 1
    stag_mg_rel_tol = 1.d-9

    gmres_rel_tol = 1.d-9
    gmres_abs_tol = 0.d0
    gmres_verbose = 1
    gmres_max_outer = 20
    gmres_max_inner = 5
    gmres_max_iter = 100
    gmres_min_iter = 1

    need_inputs = .true.

    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_gmres)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! stuff that can be read in from the command line by appending, e.g., "--precon_type 1"
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--precon_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) precon_type

       case ('--visc_schur_approx')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) visc_schur_approx

       case ('--p_norm_weight')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) p_norm_weight

       case ('--scale_factor')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) scale_factor

       case ('--mg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_verbose

       case ('--cg_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) cg_verbose

       case ('--mg_max_vcycles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_max_vcycles

       case ('--mg_minwidth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_minwidth

       case ('--mg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_bottom_solver

       case ('--mg_nsmooths_down')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_down

       case ('--mg_nsmooths_up')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_up

       case ('--mg_nsmooths_bottom')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_nsmooths_bottom

       case ('--mg_max_bottom_nlevels')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_max_bottom_nlevels

       case ('--mg_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) mg_rel_tol

       case ('--stag_mg_verbosity')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_verbosity

       case ('--stag_mg_max_vcycles')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_max_vcycles

       case ('--stag_mg_minwidth')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_minwidth

       case ('--stag_mg_bottom_solver')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_bottom_solver

       case ('--stag_mg_nsmooths_down')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_down

       case ('--stag_mg_nsmooths_up')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_up

       case ('--stag_mg_nsmooths_bottom')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_nsmooths_bottom

       case ('--stag_mg_max_bottom_nlevels')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_max_bottom_nlevels

       case ('--stag_mg_omega')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_omega

       case ('--stag_mg_smoother')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) stag_mg_smoother

       case ('--gmres_rel_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_rel_tol

       case ('--gmres_abs_tol')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_rel_tol

       case ('--gmres_verbose')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_verbose

       case ('--gmres_max_outer')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_outer

       case ('--gmres_max_inner')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_inner
      
       case ('--gmres_max_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_max_iter

      case ('--gmres_min_iter')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) gmres_min_iter

       case ('--')
          farg = farg + 1
          exit

       case default

       end select

       farg = farg + 1
    end do
    
  end subroutine probin_gmres_init

end module probin_gmres_module

