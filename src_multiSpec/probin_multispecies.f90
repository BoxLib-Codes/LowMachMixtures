module probin_multispecies_module

  use bl_types
  use bl_space
 
  implicit none

  integer, parameter :: max_species=10
  integer, parameter :: max_element=max_species*(max_species-1)/2  

  integer, save      :: nspecies,inverse_type,timeinteg_type,temp_type,chi_iterations
  real(kind=dp_t)    :: start_time
  real(kind=dp_t)    :: T_init(2) 
  real(kind=dp_t)    :: Dbar(max_element)
  real(kind=dp_t)    :: Dtherm(max_species)
  real(kind=dp_t)    :: H_offdiag(max_element)
  real(kind=dp_t)    :: H_diag(max_species)
  real(kind=dp_t)    :: fraction_tolerance
  logical            :: correct_flux,print_error_norms,plot_stag
  logical            :: is_nonisothermal,is_ideal_mixture,use_lapack
  real(kind=dp_t)    :: c_init(2,max_species)
  real(kind=dp_t)    :: c_bc(3,2,max_species)
  real(kind=dp_t)    :: alpha1,beta,delta,sigma     ! manufactured solution parameters populated in init
  
  namelist /probin_multispecies/ nspecies
  namelist /probin_multispecies/ fraction_tolerance ! For roundoff errors in mass and mole fractions
  namelist /probin_multispecies/ start_time
  namelist /probin_multispecies/ inverse_type       ! Only for LAPACK:  1=inverse, 2=pseudo inverse
  namelist /probin_multispecies/ timeinteg_type   
  namelist /probin_multispecies/ correct_flux       ! Manually ensure mass is conserved to roundoff 
  namelist /probin_multispecies/ print_error_norms   
  namelist /probin_multispecies/ is_ideal_mixture   ! If T assume Gamma=I (H=0) and simplify
  namelist /probin_multispecies/ is_nonisothermal   ! If T Soret effect will be included
  namelist /probin_multispecies/ use_lapack         ! Use LAPACK or iterative method for diffusion matrix (recommend False)
  namelist /probin_multispecies/ chi_iterations     ! number of iterations used in Dbar2chi_iterative
  namelist /probin_multispecies/ T_init     ! initial values for temperature (bottom/top, inside/outside circle, etc.)
  namelist /probin_multispecies/ temp_type  ! for initializing temperature
  namelist /probin_multispecies/ c_init   ! initial values for c
  namelist /probin_multispecies/ c_bc     ! c_i boundary conditions (dir,lohi,species)
  ! These are lower-triangules of symmetric matrices represented as vectors
  ! Number of elements is (nspecies*(nspecies-1)/2)
  ! The values are red row by row starting from top going down (this allows easy addition/deletion of new species/rows)
  ! So D_12; D_13, D_23; D_14, D_24, D_34; ...
  namelist /probin_multispecies/ Dbar       ! Maxwell-Stefan diffusion constant  
  namelist /probin_multispecies/ Dtherm     ! thermo-diffusion coefficients, only differences among elements matter
  namelist /probin_multispecies/ H_offdiag
  namelist /probin_multispecies/ H_diag     ! Diagonal of H=d^2F/dx^2, these are vectors of length nspecies
  namelist /probin_multispecies/ plot_stag  ! plot staggered velocities in separate plotfile

contains

  subroutine probin_multispecies_init()

    use f2kcli
    use parallel
    use bl_IO_module
    use bl_prof_module
    use bl_error_module
    use bl_constants_module
    use cluster_module
    
    integer            :: narg, farg
    character(len=128) :: fname
    integer            :: un
    logical            :: lexist,need_inputs
    
    narg = command_argument_count()

    ! here we set some random values to be replaced from the input file
    nspecies           = 2 
    fraction_tolerance = 1e-13 
    start_time         = 0.0d0 
    inverse_type       = 1
    timeinteg_type     = 1
    correct_flux       = .true.
    print_error_norms  = .true.
    is_ideal_mixture   = .true.
    is_nonisothermal   = .true.
    use_lapack         = .false.
    chi_iterations     = 10
    T_init             = 1.0d0
    temp_type          = 0
    c_init             = 1.0d0
    c_bc               = 0.d0
    Dbar               = 1.0d0
    Dtherm             = 0.0d0
    H_offdiag          = 0.0d0
    H_diag             = 0.0d0
    plot_stag          = .false.
 
    ! read from input file 
    need_inputs = .true.
    farg = 1
    if ( need_inputs .AND. narg >= 1 ) then
       call get_command_argument(farg, value = fname)
       inquire(file = fname, exist = lexist )
       if ( lexist ) then
          farg = farg + 1
          un = unit_new()
          open(unit=un, file = fname, status = 'old', action = 'read')
          read(unit=un, nml = probin_multispecies)
          close(unit=un)
          need_inputs = .false.
       end if
    end if

    ! also can be read in from the command line by appending 
    do while ( farg <= narg )
       call get_command_argument(farg, value = fname)
       select case (fname)

       case ('--nspecies')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) nspecies

       case ('--fraction_tolerance')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) fraction_tolerance

       case ('--start_time')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) start_time

       case ('--inverse_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) inverse_type

       case ('--timeinteg_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) timeinteg_type

       case ('--correct_flux')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) correct_flux

       case ('--print_error_norms')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) print_error_norms

       case ('--is_ideal_mixture')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) is_ideal_mixture

       case ('--is_nonisothermal')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) is_nonisothermal

       case ('--use_lapack')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) use_lapack

       case ('--chi_iterations')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) chi_iterations

       case ('--temp_type')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) temp_type

       case ('--plot_stag')
          farg = farg + 1
          call get_command_argument(farg, value = fname)
          read(fname, *) plot_stag

       case ('--')
          farg = farg + 1
          exit
       case default

       end select
       farg = farg + 1
    end do
    
    ! check that nspecies<=max_species, otherwise abort with error message
    if(nspecies.gt.max_species) then 
       call bl_error(" nspecies greater than max_species - Aborting")
       stop
    end if
    
  end subroutine probin_multispecies_init

end module probin_multispecies_module
