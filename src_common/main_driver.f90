subroutine main_driver()

  use boxlib
  use bl_IO_module
  use ml_layout_module
  use init_lowmach_module
  use init_temp_module
  use compute_mixture_properties_module
  use initial_projection_module
  use write_plotfileLM_module
  use advance_timestep_overdamped_module
  use advance_timestep_inertial_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use analysis_module
  use div_and_grad_module
  use eos_check_module
  use estdt_module
  use stag_mg_layout_module
  use macproject_module
  use stochastic_mass_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use fill_umac_ghost_cells_module
  use fill_rho_ghost_cells_module
  use ParallelRNGs 
  use mass_flux_utilities_module
  use compute_HSE_pres_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use convert_m_to_umac_module
  use sum_momenta_module
  use restart_module
  use checkpoint_module
  use project_onto_eos_module
  use probin_common_module, only: prob_lo, prob_hi, n_cells, dim_in, max_grid_size, &
                                  plot_int, chk_int, seed, stats_int, bc_lo, bc_hi, restart, &
                                  probin_common_init, print_int, project_eos_int, &
                                  advection_type, fixed_dt, max_step, cfl, &
                                  algorithm_type, variance_coef_mom, initial_variance, &
                                  barodiffusion_type
  use probin_multispecies_module, only: nspecies, Dbar, &
                                        start_time, &
                                        probin_multispecies_init
  use probin_gmres_module, only: probin_gmres_init

  use fabio_module

  implicit none

  ! quantities will be allocated with dm components
  integer, allocatable :: lo(:), hi(:)

  ! quantities will be allocated with (nlevs,dm) components
  real(kind=dp_t), allocatable :: dx(:,:)
  real(kind=dp_t)              :: dt,time,runtime1,runtime2,Dbar_max,dt_diffusive
  integer                      :: n,nlevs,i,dm,istep,ng_s,init_step,n_Dbar
  type(box)                    :: bx
  type(ml_boxarray)            :: mba
  type(ml_layout)              :: mla
  type(bc_tower)               :: the_bc_tower
  logical, allocatable         :: pmask(:)
  
  ! will be allocated on nlevels
  type(multifab), allocatable  :: rho_old(:)
  type(multifab), allocatable  :: rhotot_old(:)
  type(multifab), allocatable  :: rho_new(:)
  type(multifab), allocatable  :: rhotot_new(:)
  type(multifab), allocatable  :: Temp(:)
  type(multifab), allocatable  :: Temp_ed(:,:)
  type(multifab), allocatable  :: diff_mass_fluxdiv(:)
  type(multifab), allocatable  :: stoch_mass_fluxdiv(:)
  type(multifab), allocatable  :: umac(:,:)
  type(multifab), allocatable  :: mtemp(:,:)
  type(multifab), allocatable  :: rhotot_fc(:,:)
  type(multifab), allocatable  :: gradp_baro(:,:)
  type(multifab), allocatable  :: pi(:)
  type(multifab), allocatable  :: eta(:)
  type(multifab), allocatable  :: eta_ed(:,:)
  type(multifab), allocatable  :: kappa(:)
  type(multifab), allocatable  :: conc(:)

  !==============================================================
  ! Initialization
  !==============================================================

  call probin_common_init()
  call probin_multispecies_init() 
  call probin_gmres_init()
  
  ! Initialize random numbers *after* the global (root) seed has been set:
  call SeedParallelRNG(seed)

  ! in this example we fix nlevs to be 1
  ! for adaptive simulations where the grids change, cells at finer
  ! resolution don't necessarily exist depending on your tagging criteria,
  ! so max_levs isn't necessary equal to nlevs
  nlevs = 1

  ! dimensionality is set in inputs file
  dm = dim_in
 
  ! now that we have dm, we can allocate these
  allocate(lo(dm),hi(dm))
  allocate(rho_old(nlevs),rhotot_old(nlevs),pi(nlevs))
  allocate(rho_new(nlevs),rhotot_new(nlevs))
  allocate(Temp(nlevs),diff_mass_fluxdiv(nlevs),stoch_mass_fluxdiv(nlevs))
  allocate(umac(nlevs,dm),mtemp(nlevs,dm),rhotot_fc(nlevs,dm),gradp_baro(nlevs,dm))
  allocate(eta(nlevs),kappa(nlevs),conc(nlevs))
  if (dm .eq. 2) then
     allocate(eta_ed(nlevs,1))
     allocate(Temp_ed(nlevs,1))
  else if (dm .eq. 3) then
     allocate(eta_ed(nlevs,3))
     allocate(Temp_ed(nlevs,3))
  end if

  ! set grid spacing at each level
  ! the grid spacing is the same in each direction
  allocate(dx(nlevs,dm))
  dx(1,1:dm) = (prob_hi(1:dm)-prob_lo(1:dm)) / n_cells(1:dm)
  select case (dm) 
  case(2)
     if (dx(1,1) .ne. dx(1,2)) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy')
     end if
  case(3)
     if ((dx(1,1) .ne. dx(1,2)) .or. (dx(1,1) .ne. dx(1,3))) then
        call bl_error('ERROR: main_driver.f90, we only support dx=dy=dz')
     end if
  case default
     call bl_error('ERROR: main_driver.f90, dimension should be only equal to 2 or 3')
  end select

  ! use refined dx for next level
  do n=2,nlevs
     dx(n,:) = dx(n-1,:) / mba%rr(n-1,:)
  end do

  ! build pmask
  allocate(pmask(dm))
  pmask = .false.
  do i=1,dm
     if (bc_lo(i) .eq. PERIODIC .and. bc_hi(i) .eq. PERIODIC) then
        pmask(i) = .true.
     end if
  end do

  if (advection_type .eq. 0) then
     ng_s = 2 ! centered advection
  else if (advection_type .le. 3) then
     ng_s = 3 ! bilinear bds or unlimited quadratic bds
  else if (advection_type .eq. 4) then
     ng_s = 4 ! limited quadratic bds
  end if

  if (restart .ge. 0) then

     init_step = restart + 1

     ! build the ml_layout
     ! read in time and dt from checkpoint
     ! build and fill rho, rhotot, pi, and umac
     call initialize_from_restart(mla,time,dt,rho_old,rhotot_old,pi, &
                                  diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                  umac,pmask)

  else

     init_step = 1
     time = start_time
     
     ! tell mba how many levels and dimensionality of problem
     call ml_boxarray_build_n(mba,nlevs,dm)

     ! tell mba about the ref_ratio between levels
     ! mba%rr(n-1,i) is the refinement ratio between levels n-1 and n in direction i
     ! we use refinement ratio of 2 in every direction between all levels
     do n=2,nlevs
        mba%rr(n-1,:) = 2
     end do

     ! create a box from (0,0) to (n_cells-1,n_cells-1)
     lo(1:dm) = 0
     hi(1:dm) = n_cells(1:dm)-1
     bx = make_box(lo,hi)

     ! tell mba about the problem domain at every level
     mba%pd(1) = bx
     do n=2,nlevs
        mba%pd(n) = refine(mba%pd(n-1),mba%rr((n-1),:))
     end do

     ! initialize the boxarray at level 1 to be one single box
     call boxarray_build_bx(mba%bas(1),bx)

     ! overwrite the boxarray at level 1 to respect max_grid_size
     call boxarray_maxsize(mba%bas(1),max_grid_size)

     ! now build the boxarray at other levels
     if (nlevs .ge. 2) then
        call bl_error("Need to build boxarray for n>1")
     end if

     ! build the ml_layout, mla
     call ml_layout_build(mla,mba,pmask)

     ! don't need this anymore - free up memory
     call destroy(mba)

     do n=1,nlevs
        call multifab_build(rho_old(n)   ,mla%la(n),nspecies,ng_s)
        call multifab_build(rhotot_old(n),mla%la(n),1       ,ng_s)
        ! pi - need 1 ghost cell since we calculate its gradient
        call multifab_build(pi(n)      ,mla%la(n),1       ,1)
        call multifab_build(diff_mass_fluxdiv(n), mla%la(n),nspecies,0) 
        call multifab_build(stoch_mass_fluxdiv(n),mla%la(n),nspecies,0) 
        do i=1,dm
           call multifab_build_edge(umac(n,i),mla%la(n),1,1,i)
        end do
     end do

  end if

  deallocate(pmask)

  !=======================================================
  ! Setup boundary condition bc_tower
  !=======================================================
 
  ! bc_tower structure in memory
  ! 1:dm = velocity
  ! dm+1 = pressure
  ! dm+2 = scal_bc_comp = rhotot
  ! scal_bc_comp+1 = c_i
  ! scal_bc_comp+nspecies+1 = molfrac or massfrac (dimensionless fractions)
  ! scal_bc_comp+2*nspecies+1 = temp_bc_comp = temperature
  ! scal_bc_comp+2*nspecies+2 = tran_bc_comp = diffusion coefficients (eta,kappa,chi)
  ! It may be better if each transport coefficient has its own BC code?
  ! I think the only place this is used is average_cc_to_node/face/edge
  ! I cannot right now foresee a case where different values would be used in different places
  ! so it is OK to keep num_tran_bc_in=1. But note the same code applies to eta,kappa and chi's
  call initialize_bc(the_bc_tower,nlevs,dm,mla%pmask, &
                     num_scal_bc_in=2*nspecies+2, &
                     num_tran_bc_in=1)

  do n=1,nlevs
     ! define level n of the_bc_tower
     call bc_tower_level_build(the_bc_tower,n,mla%la(n))
  end do

  ! these quantities are populated here and defined in bc.f90
  c_bc_comp   = scal_bc_comp + 1
  mol_frac_bc_comp   = scal_bc_comp + nspecies + 1
  temp_bc_comp       = scal_bc_comp + 2*nspecies + 1

  ! build layouts for staggered multigrid solver and macproject within preconditioner
  call stag_mg_layout_build(mla)
  call mgt_macproj_precon_build(mla,dx,the_bc_tower)

  if (restart .lt. 0) then

     ! initialize rho
     call init_rho_and_umac(mla,rho_old,umac,dx,time,the_bc_tower%bc_tower_array)

     ! initialize pressure
     do n=1,nlevs
        call multifab_setval(pi(n),0.d0,all=.true.)
     end do

  end if

  do n=1,nlevs
     call multifab_build(conc(n),mla%la(n),nspecies,ng_s)
  end do

  ! compute rhotot from rho in VALID REGION
  call compute_rhotot(mla,rho_old,rhotot_old)

  ! rho to conc - NO GHOST CELLS
  call convert_rhoc_to_c(mla,rho_old,rhotot_old,conc,.true.)
  call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

  ! fill ghost cells
  do n=1,nlevs
     ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
     call multifab_fill_boundary(pi(n))
     ! fill non-periodic domain boundary ghost cells
     call multifab_physbc(pi(n),1,pres_bc_comp,1, &
                          the_bc_tower%bc_tower_array(n),dx_in=dx(n,:))
  end do

  do n=1,nlevs
     call fill_rho_ghost_cells(conc(n),rhotot_old(n),the_bc_tower%bc_tower_array(n))
  end do

  ! conc to rho - INCLUDING GHOST CELLS
  call convert_rhoc_to_c(mla,rho_old,rhotot_old,conc,.false.)

  do n=1,nlevs
     call multifab_destroy(conc(n))
  end do

  do n=1,nlevs
     do i=1,dm
        call multifab_build_edge     (mtemp(n,i),mla%la(n),1,0,i)
        call multifab_build_edge( rhotot_fc(n,i),mla%la(n),1,0,i)
        call multifab_build_edge(gradp_baro(n,i),mla%la(n),1,0,i)
     end do
  end do

  if (print_int .gt. 0) then
     if (parallel_IOProcessor()) write(*,*) "Initial state:"  
     call sum_mass(rho_old, 0) ! print out the total mass to check conservation
     ! compute rhotot on faces
     call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                             the_bc_tower%bc_tower_array)
     ! compute momentum
     call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
     call sum_momenta(mla,mtemp)
     call eos_check(mla,rho_old)
  end if

  !=======================================================
  ! Build multifabs for all the variables
  !=======================================================

  ! build multifab with nspecies component and one ghost cell
  do n=1,nlevs 
     call multifab_build(rho_new(n),           mla%la(n),nspecies,ng_s)
     call multifab_build(rhotot_new(n),        mla%la(n),1,       ng_s) 
     call multifab_build(Temp(n),              mla%la(n),1,       ng_s)
     call multifab_build(eta(n)  ,mla%la(n),1,1)
     call multifab_build(kappa(n),mla%la(n),1,1)

     ! eta and Temp on nodes (2d) or edges (3d)
     if (dm .eq. 2) then
        call multifab_build_nodal(eta_ed(n,1),mla%la(n),1,0)
        call multifab_build_nodal(Temp_ed(n,1),mla%la(n),1,0)
     else
        nodal_temp(1) = .true.
        nodal_temp(2) = .true.
        nodal_temp(3) = .false.
        call multifab_build(eta_ed(n,1),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,1),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .true.
        nodal_temp(2) = .false.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,2),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,2),mla%la(n),1,0,nodal_temp)
        nodal_temp(1) = .false.
        nodal_temp(2) = .true.
        nodal_temp(3) = .true.
        call multifab_build(eta_ed(n,3),mla%la(n),1,0,nodal_temp)
        call multifab_build(Temp_ed(n,3),mla%la(n),1,0,nodal_temp)
     end if

  end do

  ! allocate and build multifabs that will contain random numbers
  if (algorithm_type .eq. 0 .or. algorithm_type .eq. 1) then
     n_rngs = 1
  else if (algorithm_type .eq. 2) then
     n_rngs = 2
  end if
  call init_mass_stochastic(mla,n_rngs)
  call init_m_stochastic(mla,n_rngs)

  ! fill random flux multifabs with new random numbers
  call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)
  call fill_m_stochastic(mla)

  !=====================================================================
  ! Initialize values
  !=====================================================================

  ! initialize Temp
  call init_Temp(Temp,dx,time,the_bc_tower%bc_tower_array)
  if (dm .eq. 2) then
     call average_cc_to_node(nlevs,Temp,Temp_ed(:,1),1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  else if (dm .eq. 3) then
     call average_cc_to_edge(nlevs,Temp,Temp_ed,1,tran_bc_comp,1,the_bc_tower%bc_tower_array)
  end if

  if (barodiffusion_type .gt. 0) then

     ! this computes an initial guess at p using HSE
     call compute_HSE_pres(mla,rhotot_old,pi,dx,the_bc_tower)

     ! compute grad p for barodiffusion
     call compute_grad(mla,pi,gradp_baro,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

  end if

  ! initialize eta and kappa
  call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_old,rhotot_old,Temp,dx, &
                         the_bc_tower%bc_tower_array)

  call fill_umac_ghost_cells(mla,umac,eta_ed,dx,time,the_bc_tower)

  if (restart .lt. 0) then

     ! add initial momentum fluctuations - only call in inertial code for now
     ! Note, for overdamped code, the steady Stokes solver will wipe out the initial condition
     if (algorithm_type .eq. 0 .and. initial_variance .ne. 0.d0) then
        call add_m_fluctuations(mla,dx,initial_variance*variance_coef_mom, &
                                umac,rhotot_old,Temp,the_bc_tower)
     end if

     if (fixed_dt .gt. 0.d0) then
        dt = fixed_dt
     else
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        dt = min(dt,dt_diffusive)
     end if
     
  end if

  !=====================================================================
  ! Process initial conditions
  !=====================================================================

  if (restart .lt. 0) then
     
     ! initial projection - only truly needed for inertial algorithm
     ! for the overdamped algorithm, this only changes the reference state for the first
     ! gmres solve in the first time step
     ! Yes, I think in the purely overdamped version this can be removed
     ! In either case the first ever solve cannot have a good reference state
     ! so in general there is the danger it will be less accurate than subsequent solves
     ! but I do not see how one can avoid that
     ! From this perspective it may be useful to keep initial_projection even in overdamped
     ! because different gmres tolerances may be needed in the first step than in the rest
     if (algorithm_type .eq. 0) then
        call initial_projection(mla,umac,rho_old,rhotot_old,gradp_baro,diff_mass_fluxdiv, &
                                stoch_mass_fluxdiv, &
                                Temp,eta,eta_ed,dt,dx,the_bc_tower)
     end if

     if (print_int .gt. 0) then
        if (parallel_IOProcessor()) write(*,*) "After initial projection:"  
        call sum_mass(rho_old,0) ! print out the total mass to check conservation
        ! compute rhotot on faces
        call average_cc_to_face(nlevs,rhotot_old,rhotot_fc,1,scal_bc_comp,1, &
                                the_bc_tower%bc_tower_array)
        ! compute momentum
        call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
        call sum_momenta(mla,mtemp)
        call eos_check(mla,rho_old)
     end if   

     ! write initial plotfile
     if (plot_int .gt. 0) then
        if (parallel_IOProcessor()) then
           write(*,*), 'writing initial plotfile 0'
        end if
        call write_plotfileLM(mla,"plt",rho_old,rhotot_old,Temp,umac,pi,0,dx,time)
     end if
     
     ! print out projection (average) and variance)
     if (stats_int .gt. 0) then
        call print_stats(mla,dx,0,time,umac=umac,rho=rho_old,temperature=Temp)
     end if

  end if

  !=======================================================
  ! Begin time stepping loop
  !=======================================================

  do istep=init_step,max_step

     if (fixed_dt .le. 0.d0) then
        call estdt(mla,umac,dx,dt)
        n_Dbar = nspecies*(nspecies-1)/2
        Dbar_max = maxval(Dbar(1:n_Dbar))
        dt_diffusive = cfl*dx(1,1)**2/(2*dm*Dbar_max)
        dt = min(dt,dt_diffusive)
     end if

     if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"Begin Advance; istep =",istep,"dt =",dt,"time =",time
     end if

     runtime1 = parallel_wtime()

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! advance the solution by dt
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! notes: eta, eta_ed, and kappa could be built and initialized within the advance routines
      ! but for now we pass them around (it does save a few flops)
      ! diff/stoch_mass_fluxdiv could be built locally within the overdamped
      ! routine, but since we have them around anyway for inertial we pass them in
      if (algorithm_type .eq. 0) then
         call advance_timestep_inertial(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                        gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                        diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                        dx,dt,time,the_bc_tower,istep)
      else if (algorithm_type .eq. 1 .or. algorithm_type .eq. 2) then
         call advance_timestep_overdamped(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                          gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                          diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                          dx,dt,time,the_bc_tower,istep)
      end if

      time = time + dt

      if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,"End Advance; istep =",istep,"DT =",dt,"TIME =",time
      end if

      runtime2 = parallel_wtime() - runtime1
      call parallel_reduce(runtime1, runtime2, MPI_MAX, proc=parallel_IOProcessorNode())
      if (parallel_IOProcessor()) then
        if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) ) &
           print*,'Time to advance timestep: ',runtime1,' seconds'
      end if
      
      if ( (print_int .gt. 0 .and. mod(istep,print_int) .eq. 0) &
          .or. &
          (istep .eq. max_step) ) then
          if (parallel_IOProcessor()) write(*,*) "At time step ", istep, " t=", time           
          call sum_mass(rho_new, istep) ! print out the total mass to check conservation
          ! compute rhotot on faces
          call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,scal_bc_comp,1, &
                                  the_bc_tower%bc_tower_array)
          ! compute momentum
          call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)
          call sum_momenta(mla,mtemp)
          call eos_check(mla,rho_new)
      end if

     ! project rho and rho1 back onto EOS
     if ( project_eos_int .gt. 0 .and. mod(istep,project_eos_int) .eq. 0) then
        call project_onto_eos(mla,rho_new)
     end if

      ! We do the analysis first so we include the initial condition in the files if n_steps_skip=0
      if (istep >= n_steps_skip) then

         ! write plotfile at specific intervals
         if (plot_int.gt.0 .and. ( (mod(istep,plot_int).eq.0) .or. (istep.eq.max_step)) ) then
            if (parallel_IOProcessor()) then
               write(*,*), 'writing plotfiles at timestep =', istep 
            end if
            call write_plotfileLM(mla,"plt",rho_new,rhotot_new,Temp,umac,pi,istep,dx,time)
         end if

         ! write checkpoint at specific intervals
         if ((chk_int.gt.0 .and. mod(istep,chk_int).eq.0)) then
            if (parallel_IOProcessor()) then
               write(*,*), 'writing checkpoint at timestep =', istep 
            end if
            call checkpoint_write(mla,rho_new,rhotot_new,pi,diff_mass_fluxdiv, &
                                  stoch_mass_fluxdiv,umac,time,dt,istep)
         end if

         ! print out projection (average) and variance
!         if ( (stats_int > 0) .and. &
!               (mod(istep,stats_int) .eq. 0) ) then
!            ! Compute vertical and horizontal averages (hstat and vstat files)   
!            call print_stats(mla,dx,istep,time,umac=umac,rho=rho_new,temperature=Temp)            
!         end if

      end if

      ! set old state to new state
      do n=1,nlevs
         call multifab_copy_c(rho_old(n)   ,1,rho_new(n)   ,1,nspecies,rho_old(n)%ng)
         call multifab_copy_c(rhotot_old(n),1,rhotot_new(n),1,1       ,rhotot_old(n)%ng)
      end do

  end do

  !=======================================================
  ! Destroy multifabs and layouts
  !=======================================================

  call destroy_mass_stochastic(mla)
  call destroy_m_stochastic(mla)

  do n=1,nlevs
     call multifab_destroy(rho_old(n))
     call multifab_destroy(rhotot_old(n))
     call multifab_destroy(rho_new(n))
     call multifab_destroy(rhotot_new(n))
     call multifab_destroy(Temp(n))
     call multifab_destroy(diff_mass_fluxdiv(n))
     call multifab_destroy(stoch_mass_fluxdiv(n))
     call multifab_destroy(pi(n))
     call multifab_destroy(eta(n))
     call multifab_destroy(kappa(n))
     do i=1,dm
        call multifab_destroy(umac(n,i))
        call multifab_destroy(mtemp(n,i))
        call multifab_destroy(rhotot_fc(n,i))
        call multifab_destroy(gradp_baro(n,i))
     end do
     do i=1,size(eta_ed,dim=2)
        call multifab_destroy(eta_ed(n,i))
        call multifab_destroy(Temp_ed(n,i))
     end do
  end do
  deallocate(lo,hi,dx)
  deallocate(rho_old,rhotot_old,Temp,umac)
  deallocate(diff_mass_fluxdiv,stoch_mass_fluxdiv)
  call stag_mg_layout_destroy()
  call mgt_macproj_precon_destroy()
  call destroy(mla)
  call bc_tower_destroy(the_bc_tower)

end subroutine main_driver
