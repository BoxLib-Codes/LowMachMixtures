!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! HydroGrid analysis for spectra/rendering and other serial routines
  ! Make a copy of all the data on one of the processors and then call the serial routines
  ! Important note: RESTARTS NOT IMPLEMENTED in analysis code!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module analyze_spectra_module

  use fabio_module
  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use HydroGridModule
  use HydroGridCInterface 
  use probin_common_module, only: n_cells, prob_lo, prob_hi, &
       hydro_grid_int, project_dir, max_grid_projection, stats_int, n_steps_save_stats, &
       center_snapshots, analyze_conserved

  implicit none

  private
  public :: initialize_hydro_grid, analyze_hydro_grid, save_hydro_grid, finalize_hydro_grid, &
       print_stats

  ! For statistical analysis of fluctuating fields
  type (HydroGrid), target, save :: grid, grid_2D, grid_1D
  
  ! We collect all the data on processor 1 for the analysis stuff due to the FFTs etc:
  integer, save :: ncells(3), ncells2D(3), ncells1D(3)
  ! These are ordered so that velocities come first, then densities, then temperature, and lastly scalars
  ! Any of these can be omitted (not present)
  type(multifab), save ::  s_serial, s_dir, s_projected, s_var ! Only one level here
  ! Must keep the layout around until the end!
  type(layout), save :: la_serial, la_dir, la_projected

  ! number of (non-velocity) scalars to be analyzed
  integer, save :: nvar=1, nscal_analysis=1, nspecies_analysis=1
      ! nvar is total number of variables (vector and scalar)
      ! nscal_analysis is total number of scalar variables to analyze
  logical, save :: exclude_last_species=.true. ! Whether to use format 
      ! true=(rho,rho_1,...,rho_{n-1}) or false=(rho_1,...,rho_{n})

  ! Note: analyze_conserved determines whether we analyze densities rho_k or mass fractions w_k=rho_k/rho
  ! All other quantities are assumed to be primitive
  ! but of course you can pass conserved ones and the code won't know you lied!
  
  ! project_dir determines a direction for performing projections, if desired
  ! If positive, it means only analyze projections (in parallel, much more efficient)
  ! If negative, it does serialized analysis of the whole grid *and* the projections
  ! If zero, it only does serialized analysis of the whole grid
   
contains   

  ! When we call HydroGrid we tall it there an nspecies+1 species so that it analyzes all species instead of omitting the last one
  ! This is redundant but conventent if one is also interested in the last species
  ! Pass nspecies_in=0 if there are no compositional variables
  subroutine initialize_hydro_grid(mla,s_in,dt,dx,namelist_file, &
                                   nspecies_in, nscal_in, exclude_last_species_in, &
                                   analyze_velocity, analyze_density, analyze_temperature, &
                                   heat_capacity_in)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_in(:) ! A cell-centered multifab on the desired grid (to grab layout and grid size from)
    real(dp_t)     , intent(inout) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: namelist_file ! Where to read the namelists from
    integer        , intent(in   ) :: nspecies_in ! Number of species (number of partial densities)
       ! Note: Pass nspecies_in=0 if there are no compositional variables
    integer        , intent(in   ) :: nscal_in
      ! Additional number of scalar fields if any
    logical, intent(in) :: exclude_last_species_in, analyze_velocity, analyze_density, analyze_temperature
        ! Pass exclude_last_species=.false. if you want to analyze all nspecies densities/concentrations
        ! Or pass exclude_last_species=.true. if you want to pass total density rho as the first scalar and 
        !    then only include only n_species-1 additional partial densities/concentrations
    real(dp_t), intent(in), optional :: heat_capacity_in(nspecies_in) ! This assumes constant heat capacity (ideal gas mixture)   

    ! local
    type(box)  :: bx_serial, bx_dir, bx_projected
    type(boxarray)  :: bxa_serial, bxa_dir, bxa_projected
    integer, dimension(mla%dim) :: lo, hi
    integer :: ncells_tmp(nMaxDims)
    integer :: nlevs, dm, pdim, max_grid_dir(3), max_grid_projected(3)
    ! The analysis codes always works in 3D
    real(dp_t) :: grid_dx(3)
    real(dp_t) :: heat_capacity(nspecies_in+1)

    nlevs = mla%nlevel
    dm = mla%dim

    nspecies_analysis = nspecies_in ! Number of concentration variables
    exclude_last_species = exclude_last_species_in
    
    if(present(heat_capacity_in)) then
      heat_capacity(1:nspecies_in) = heat_capacity_in
      heat_capacity(nspecies_in+1) = 0.0d0
    else
      heat_capacity = 1.0_dp_t ! Default value
    end if

    if(nlevs>1) call parallel_abort("HydroGrid analysis only implemented for a single level!")

    ! Calculate total number of scalar variables being analyzed to allocate storage
    nvar=nscal_in ! nvar holds the total number of variables to analyze (scalar and vector)
    if(analyze_density)  nvar = nvar + nspecies_analysis + 1 ! We always store total density
    if(analyze_temperature) nvar = nvar + 1
    nscal_analysis = nvar ! Total number of scalar variables to analyze
    if(analyze_velocity) nvar = nvar + dm
    if ( parallel_IOprocessor() ) then
       write(*,*) "Allocating storage for ", nvar, " variables to analyze using HydroGrid"
    end if

    ncells(1) = n_cells(1)
    ncells(2) = n_cells(2)
    ncells(3) = 1
    if(dm>2) ncells(3) = n_cells(3)

    pdim=abs(project_dir) ! Axes to project along
    if(pdim>0) then
       ncells2D=ncells
       ncells2D(pdim)=1   
       ncells1D=1
       ncells1D(pdim)=ncells(pdim)
    end if

    if ((stats_int > 0) .or. (project_dir > 0)) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! build s_dir (multifab with "tall skinny boxes")
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! It is hard to write a general algorithm for figuring out the best way to split into boxes
       ! So we make the user specify that via max_grid_projection
       select case(pdim)
       case(1)
          max_grid_dir=(/n_cells(1),max_grid_projection(1),max_grid_projection(2)/)
       case(2)
          max_grid_dir=(/max_grid_projection(1),n_cells(2),max_grid_projection(2)/)
       case(3)
          max_grid_dir=(/max_grid_projection(1),max_grid_projection(2),n_cells(3)/)
       case default
          call bl_error("project_dir must be between 1 and 3 if stats_int>0")
       end select

       bx_dir = s_in(1)%la%lap%pd                   ! set bx_dir to problem domain
       call boxarray_build_bx(bxa_dir, bx_dir)      ! build a boxarray containing only one box
       call boxarray_maxsize(bxa_dir, max_grid_dir) ! chop into tall skinny boxes
       if ( parallel_IOprocessor() ) then
          call print(bxa_dir, 'Analysis tall/skinny BoxArray')
          print*,''
       end if
       call layout_build_ba (la_dir, bxa_dir, bx_dir, mla%pmask)
       call destroy(bxa_dir)
       call multifab_build(s_dir, la_dir, nvar, 0)

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! build s_projected (multifab with reduced dimension that holds the projection of s_dir)
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! first set lo and hi to the problem domain
       lo = 0
       hi(1) = n_cells(1)-1
       if (dm > 1) then   
          hi(2) = n_cells(2) - 1        
          if (dm > 2)  then
             hi(3) = n_cells(3) -1
          endif
       endif
       ! then reduce the dimensionality by 1
       hi(pdim)=lo(pdim)
       select case(pdim)
       case(1)
          max_grid_projected=(/1,max_grid_projection(1),max_grid_projection(2)/)
       case(2)
          max_grid_projected=(/max_grid_projection(1),1,max_grid_projection(2)/)
       case(3)
          max_grid_projected=(/max_grid_projection(1),max_grid_projection(2),1/)
       case default
          call bl_error("project_dir must be between 1 and 3")
       end select
       
       call box_build_2(bx_projected,lo,hi)                     ! set bx_projected to reduced dimension problem domain
       call boxarray_build_bx(bxa_projected,bx_projected)       ! build a boxarray containing only one box
       call boxarray_maxsize(bxa_projected, max_grid_projected) ! chop up boxes
       if ( parallel_IOprocessor() ) then
          call print(bxa_projected, 'Analysis projected BoxArray')
          print*,''
       end if
       ! force same mapping as la_dir
       call layout_build_ba(la_projected, bxa_projected, bx_projected, mla%pmask, &
                            explicit_mapping=get_proc(la_dir))
       call destroy(bxa_projected)
       call multifab_build(s_projected, la_projected, nvar, 0)

       if (stats_int > 0) then

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_var (multifab with reduced dimension that holds the variance of s_dir)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! has same layout as s_projected
          call multifab_build(s_var, la_projected, nvar, 0)

       end if

    end if

    if (abs(hydro_grid_int)>0) then

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Build serialized MultiFabs for HydroGrid analysis
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
       if (project_dir > 0) then ! Serialize projection analysis only

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_serial (serial version of s_projected)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! set bx_serial to reduced dimension problem domain
          call box_build_2(bx_serial,lo,hi)
          ! build a boxarray containing only one box
          call boxarray_build_bx(bxa_serial,bx_serial)
          if ( parallel_IOprocessor() ) then
             call print(bxa_serial, 'Analysis serial BoxArray')
             print*,''
          end if
          call layout_build_ba(la_serial, bxa_serial, bx_serial, mla%pmask, &
               explicit_mapping=(/parallel_IOProcessorNode()/) )
          call destroy(bxa_serial)
          call multifab_build(s_serial, la_serial, nvar, 0)

       else ! Serialize the whole analysis

          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! build s_serial (serial version of s_in)
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! set bx_serial to problem domain
          bx_serial = s_in(1)%la%lap%pd
          ! build a boxarray containing only one box
          call boxarray_build_bx(bxa_serial, bx_serial) 
          call layout_build_ba(la_serial, bxa_serial, bx_serial, mla%pmask, &
               explicit_mapping=(/parallel_IOProcessorNode()/) )
          if ( parallel_IOprocessor() ) then
             call print(bxa_serial, 'Analysis serial BoxArray')
             print*,''
          end if
          call destroy(bxa_serial)       
          call multifab_build(s_serial, la_serial, nvar, 0)

       end if

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Now initialize the analysis code itself:   
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       if(parallel_IOProcessor()) then
       
          grid_dx=1.0_dp_t ! Default grid spacing
          grid_dx(1:dm)=dx(1,1:dm) 

          if(project_dir <= 0) then ! Perform analysis on the full 3D grid
            ncells_tmp=nCells
          else ! Create a fake object that will not be used, with grid size 1x1x1
            ! This way the namelist input files do not have to change     
            ncells_tmp(:) = 1
          end if     
          call createHydroAnalysis (grid, nCells=ncells_tmp, nSpecies = nspecies_analysis+1, &
             isSingleFluid = .true., nVelocityDimensions = dm, nPassiveScalars = nscal_in, &
             systemLength = ncells*grid_dx, heatCapacity = heat_capacity, &
             timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
             structFactMultiplier = 1.0_dp_t )

          if(project_dir/=0) then
             ! Also perform analysis on a projected grid (averaged along project_dir axes)
             call createHydroAnalysis (grid_2D, nCells=ncells2D, nSpecies = nspecies_analysis+1, &
                  isSingleFluid = .true., nVelocityDimensions = dm, nPassiveScalars = nscal_in, &
                  systemLength = nCells*grid_dx, heatCapacity = heat_capacity, &
                  timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
                  structFactMultiplier = 1.0_dp_t )
          end if

          if(project_dir/=0) then ! Also perform analysis on a 1D grid (along project_dir only)
             call createHydroAnalysis (grid_1D, nCells=ncells1D, nSpecies = nspecies_analysis+1, &
                  isSingleFluid = .true., nVelocityDimensions = dm, nPassiveScalars = nscal_in, &
                  systemLength = nCells*grid_dx, heatCapacity = heat_capacity, &
                  timestep = abs(hydro_grid_int)*dt, fileUnit=namelist_file, &
                  structFactMultiplier = 1.0_dp_t )
          end if

       end if

    end if
      
  end subroutine

  ! Files can be numbered either by id or by step, here we use step so restarts work seemlessly
  subroutine save_hydro_grid(id, step) ! This also *resets* all counters
      integer, intent(in), optional :: id, step ! We can use either one to number files here
      
      if((.not.parallel_IOProcessor()) .or. (hydro_grid_int<=0)) return
      
      if(present(step)) write(*,*) "Saving HydroGrid analysis results at time step ", step
      
      if(project_dir<=0) call writeToFiles(grid, id=step)
      if(project_dir/=0) then
         call writeToFiles(grid_2D, id=step)
      end if
      if(project_dir/=0) then
         call writeToFiles(grid_1D, id=step)
      end if
      
  end subroutine
  
  subroutine finalize_hydro_grid()

    if ((stats_int > 0) .or. (project_dir > 0)) then
       call multifab_destroy(s_dir)
       call multifab_destroy(s_projected)
       if (stats_int > 0) then
          call multifab_destroy(s_var)
       end if
       call destroy(la_dir)
       call destroy(la_projected)
    end if
    if (abs(hydro_grid_int)>0) then
       call multifab_destroy(s_serial)
       call destroy(la_serial)
    end if
    
    if((.not.parallel_IOProcessor())) return
    
    if(hydro_grid_int<0) then ! Only do clean-up

       call destroyHydroAnalysis(grid)
       if(project_dir/=0) then
          call destroyHydroAnalysis(grid_2D) 
          call destroyHydroAnalysis(grid_1D)
       end if
    
    else if(hydro_grid_int>0) then ! Save HydroGrid data and clean-up

       if(project_dir<=0) then
         if(n_steps_save_stats <= 0) call writeToFiles(grid)
       end if  
       call destroyHydroAnalysis(grid)
       if(project_dir/=0) then
          if(n_steps_save_stats <= 0) call writeToFiles(grid_2D)
          call destroyHydroAnalysis(grid_2D) 
       end if
       if(project_dir/=0) then
          if(n_steps_save_stats <= 0) call writeToFiles(grid_1D)
          call destroyHydroAnalysis(grid_1D)
       end if
    
    end if   

  end subroutine

  subroutine analyze_hydro_grid(mla,dt,dx,step,umac,rho,temperature,scalars)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho, temperature, scalars
  
    if(project_dir>0) then ! Only do 2D analysis
      !write(*,*) "Calling analyze_hydro_grid_parallel"
      call analyze_hydro_grid_parallel(mla,dt,dx,step,umac,rho)
    else ! Analyze the whole grid
      !write(*,*) "Calling analyze_hydro_grid_serial"
      call analyze_hydro_grid_serial(mla,dt,dx,step,umac,rho,temperature,scalars)
    end if

  end subroutine

  ! This routines collects a bunch of different variables into a single multifab for easier analysis
  ! The components are ordered as, whenever present:
  ! umac, rho, rho_k (analyze_conserved=T) or c_k=rho_k/rho, T, scalars
  subroutine gather_hydro_grid(mla,s_hydro,umac,rho,temperature,scalars, variable_names)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: s_hydro ! Output: collected grid data for analysis
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho, temperature, scalars
    character(len=20), intent(out), optional :: variable_names(nvar)
  
    ! Local
    integer :: i, ii, jj, n, nlevs, dm, pdim, comp, n_passive_scals, species
    type(multifab) :: mac_cc(mla%nlevel)
 
    nlevs = mla%nlevel
    dm = mla%dim
    
    comp=1    
    if(present(umac)) then
       do n=1,nlevs
          call multifab_build(mac_cc(n), mla%la(n), dm, 0)
       end do    
       call StaggeredToCentered(mac_cc, mla, umac)
       call copy(s_hydro,1,mac_cc(1),1,dm)
       do n=1,nlevs
          call multifab_destroy(mac_cc(n))
       end do
       if(present(variable_names)) then
         variable_names(comp)="v_x"
         variable_names(comp+1)="v_y"
         if(dm>2) variable_names(comp+2)="v_z"
       end if
       comp=comp+dm
    end if
    
    ! Now all the scalars
    n_passive_scals=nscal_analysis
    if(present(rho)) then
       
       if(n_passive_scals<nspecies_analysis+1) then
         call parallel_abort("Insufficient nscal_analysis to store densities")
       end if

       if(nspecies_analysis==0) then ! Copy only rho to s_hydro
          call copy(s_hydro,comp,rho(1),1,1)
       else if(exclude_last_species) then
          call copy(s_hydro,comp,rho(1),1,nspecies_analysis) ! Copy rho,rho_1,...,rho_{n-1}
          ! Now compute rho_n as the difference rho-\sum_{i=1}^{n-1} rho_i
          call copy(s_hydro,nspecies_analysis+comp,s_hydro,comp)
          do species=1, nspecies_analysis-1
             call multifab_sub_sub_c(s_hydro,nspecies_analysis+comp,s_hydro,species+comp,1)
          end do
       else
          call copy(s_hydro,comp+1,rho(1),1,nspecies_analysis) ! copy rho_1,...,rho_n
          ! Compute total density as a sum of densities rho=\sum_{i=1}^{n} rho_i
          call setval(s_hydro, 0.0d0, comp)
          do species=1, nspecies_analysis
             call multifab_plus_plus_c(s_hydro,comp,s_hydro,species+comp,1)
          end do
       end if
       
       ! convert to mass fractions instead of densities (primitive variables)
       if(.not.analyze_conserved) then
          do species=1, nspecies_analysis
             call multifab_div_div_c(s_hydro,species+comp,s_hydro,comp,1)
          end do
       end if

       if (present(variable_names)) then
          variable_names(comp)="rho"
          if (analyze_conserved) then
             do species=1, nspecies_analysis
               write(variable_names(comp+species),"(A,I0)") "rho_", species
             end do
          else
             do species=1, nspecies_analysis
               write(variable_names(comp+species),"(A,I0)") "c_", species
             end do
          end if
       end if
    
       comp=comp+nspecies_analysis+1
       n_passive_scals=n_passive_scals-nspecies_analysis-1
    end if
    
    if(present(temperature)) then
       if(n_passive_scals<1) then
         call parallel_abort("Insufficient nscal_analysis to store temperature")
       end if
       call copy(s_hydro,comp,temperature(1),1,1)
       if(present(variable_names)) then
         variable_names(comp)="T"
       end if
       comp=comp+1    
       n_passive_scals=n_passive_scals-1
    end if
    
    if(present(scalars)) then
       if(n_passive_scals<1) then
         call parallel_abort("No room left to store passive scalars, nscal_analysis too small")
       end if
       if(present(variable_names)) then
          do species=1, n_passive_scals
            write(variable_names(comp+species-1),"(A,I0)") "s_", species
          end do
       end if
       call copy(s_hydro,comp,scalars(1),1,n_passive_scals)
    end if
  
  end subroutine

  subroutine point_to_hydro_grid(mla,variables,p_v, p_rho, p_c, p_T, p_s, &
                                 umac,rho,temperature,scalars)
    type(ml_layout), intent(in   ) :: mla
    real(kind=dp_t), dimension(:,:,:,:), target, intent(in)   :: variables
    real(kind=dp_t), pointer, dimension(:,:,:,:), intent(out) :: p_v, p_rho, p_c, p_T, p_s
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho, temperature, scalars
  
    ! Local
    integer :: i, ii, jj, n, nlevs, dm, pdim, comp, n_passive_scals, species
    
    nlevs = mla%nlevel
    dm = mla%dim

    comp=1  
    if(present(umac)) then
       p_v=>variables(:,:,:,comp:comp+dm-1)
       !write(*,*) "p_v LBOUND=", lbound(p_v), " UBOUND=", ubound(p_v)
       !write(*,*) "p_v(1,1,1,:)=", p_v(1,1,1,:)
       comp=comp+dm
    else
       p_v=>NULL()   
    end if

    ! Now all the scalars
    n_passive_scals=nscal_analysis
    ! Now gather density
    if(present(rho)) then
       p_rho => variables(:,:,:,comp:comp)
       p_c => variables(:,:,:,comp+1:comp+nspecies_analysis)
       !write(*,*) "p_rho LBOUND=", lbound(p_rho), " UBOUND=", ubound(p_rho)
       !write(*,*) "p_rho(1,1,1,1)=", p_rho(1,1,1,1)
       !write(*,*) "p_c LBOUND=", lbound(p_c), " UBOUND=", ubound(p_c)
       !write(*,*) "p_c(1,1,1,:)=", p_c(1,1,1,:)
       comp=comp+nspecies_analysis+1
       n_passive_scals=n_passive_scals-nspecies_analysis-1
    else
       p_rho=>NULL()   
       p_c=>NULL()             
    end if

    if(present(temperature)) then
       p_T => variables(:,:,:,comp:comp)
       comp=comp+1    
       n_passive_scals=n_passive_scals-1
       !write(*,*) "p_T LBOUND=", lbound(p_T), " UBOUND=", ubound(p_T)
    else
       p_T=>NULL()  
    end if

    if(present(scalars)) then
       p_s => variables(:,:,:,comp:nscal_analysis)
       !write(*,*) "p_s LBOUND=", lbound(p_s), " UBOUND=", ubound(p_s)
    else
       p_s=>NULL()   
    end if

  end subroutine
  
  ! This analysis routine is serialized:
  ! All the data from all the boxes is collected onto one processor and then analyzed
  ! The advantage of this is that one can do some customized analysis that is hard to do in parallel
  ! It is also the best one can do when one wants to analyze the whole grid instead of just projections
  subroutine analyze_hydro_grid_serial(mla,dt,dx,step,umac,rho,temperature,scalars)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho, temperature, scalars

    ! Local
    integer :: i, ii, jj, n, nlevs, dm, pdim, comp, n_passive_scals, species
    real(kind=dp_t), pointer, dimension(:,:,:,:) :: variables, p_v, p_rho, p_c, p_T, p_s
    real(kind=dp_t), dimension(:,:), allocatable :: variables_1D ! Projection of data along project_dir
    
    nlevs = mla%nlevel
    dm = mla%dim
    pdim=abs(project_dir)

    ! Gather all the hydro data into a single multifab
    call gather_hydro_grid(mla,s_serial,umac,rho,temperature,scalars)
    
    if(parallel_IOProcessor()) then  
       !write(*,*) "Calling updateHydroAnalysis on 3D+2D+1D grids (serial)"

       ! Get to the actual data
       variables  => dataptr(s_serial,1) ! Gets all of the components
       
       if(pdim>0) then

          allocate ( variables_1D(lbound(variables,dim=pdim) : ubound(variables,dim=pdim), &
            lbound(variables,dim=4) : ubound(variables,dim=4)) )

          ! Average the data in the hyperplanes perpendicular to pdim:
          do jj =  lbound(variables,dim=4), ubound(variables,dim=4)
             do ii = lbound(variables,dim=pdim), ubound(variables,dim=pdim)     
                if (pdim .eq. 1) then
                   variables_1D(ii,jj) = sum(variables(ii,:,:,jj))
                else if (pdim .eq. 2) then
                   variables_1D(ii,jj) = sum(variables(:,ii,:,jj))
                else if (pdim .eq. 3) then
                   variables_1D(ii,jj) = sum(variables(:,:,ii,jj))
                end if
             end do
          end do
          
       end if
       
       ! Point pointers to the right places, if present, NULL otherwise
       call point_to_hydro_grid(mla,variables,p_v, p_rho, p_c, p_T, p_s, &
                                umac,rho,temperature,scalars)
       if(.false.) write(*,*) "associated=", associated(p_v), associated(p_rho), &
                   associated(p_c), associated(p_T), associated(p_s)

       call updateHydroAnalysisPrimitive (grid, velocity=p_v, &
           density=p_rho, concentration=p_c, temperature=p_T, scalars=p_s)
       
       ! The following analysis we only do for concentrations for now
       if((project_dir/=0).and.(present(rho))) then ! Also average along projections
          call updateHydroAnalysisPrimitive (grid_2D, &
             density=SUM(p_rho, DIM=pdim)/grid%nCells(pdim), &
             concentration=SUM(p_c, DIM=pdim)/grid%nCells(pdim) )
          if(present(umac)) then
            comp=1+dm
          else
            comp=1
          end if       
          call updateHydroAnalysisPrimitive (grid_1D, density=variables_1D(:,comp), &
                  concentration=variables_1D(:,comp:comp+nspecies_analysis))                     
       end if

       if(pdim>0) deallocate(variables_1D)
    end if   

  end subroutine

  ! This routine first projects onto the project_dir direction and the perpendicular 
  ! hyperpane in parallel.  Then it serializes the analysis of those projections
  subroutine analyze_hydro_grid_parallel(mla,dt,dx,step,umac,rho,temperature,scalars)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in) :: dt
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho,temperature,scalars
   
    ! Local
    integer species, n_passive_scals, comp
    integer nlevs,dm,pdim,i,n,qdim,qdim1,qdim2
    integer lo(mla%dim),hi(mla%dim)
    integer ii,jj,kk

    real(kind=dp_t), pointer, dimension(:,:,:,:) :: variables, sdp, spp, svp, &
                                                    p_v, p_rho, p_c, p_T, p_s
    real(kind=dp_t), dimension(:,:,:,:), allocatable, target :: variables_1D, variables_1D_proc
      ! Projection of data along project_dir

    nlevs = mla%nlevel
    dm = mla%dim
    pdim  = abs(project_dir)

    ! qdim is orthogonal to pdim
    if (pdim .eq. 1) then
       qdim  = 2
       qdim2 = 3
    else if (pdim .eq. 2) then
       qdim  = 1
       qdim2 = 3
    else if (pdim .eq. 3) then
       qdim  = 1
       qdim2 = 2
    end if

    ! Re-distribute the full grid into a grid of "tall skinny boxes"
    ! These boxes are not distributed along project_dim so we can do local analysis easily
    ! -------------------------
    call gather_hydro_grid(mla,s_dir,umac,rho,temperature,scalars)
    
    ! Compute s_projected as the average along project_dim
    ! -------------------------
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       spp => dataptr(s_projected, i)
       lo = lwb(get_box(s_dir, i))
       hi = upb(get_box(s_dir, i))
       ! put sum ( s_dir / ncell_pdim ) into s_projected
       if (pdim .eq. 1) then
          spp(0,:,:,:)=SUM( sdp, DIM=1 ) / dble(hi(1)-lo(1)+1)
       else if (pdim .eq. 2) then
          spp(:,0,:,:)=SUM( sdp, DIM=2 ) / dble(hi(2)-lo(2)+1)
       else if (pdim .eq. 3) then
          spp(:,:,0,:)=SUM( sdp, DIM=3 ) / dble(hi(3)-lo(3)+1)
       end if
    end do

    ! Now collect the projected data from s_projected into a single box in s_serial 
    call copy(s_serial,1,s_projected,1,nvar)
    
    ! We also want to collect 1D data for the projection onto the direction pdim
    ! We need to do the averaging for each of the boxes in s_dir and then do an mpi_reduction
    allocate ( variables_1D(1:ncells(pdim), 1, 1, nvar), variables_1D_proc(1:ncells(pdim), 1, 1, nvar) )    
    call average_1D(variables_1D_proc) ! Sum over cells in each box    
    ! sum reduction: Note that dividing by the number of cells is already done in average_1D
    do n=1,nvar
       call parallel_reduce(variables_1D(:,1,1,n), variables_1D_proc(:,1,1,n), MPI_SUM, &
                            proc=parallel_IOProcessorNode())
    end do

    if(parallel_IOProcessor()) then  
       !write(*,*) "Calling updateHydroAnalysis on 2D+1D grid (parallel)"

       variables  => dataptr(s_serial,1) ! Gets all of the components

       ! Point pointers to the right places, if present, NULL otherwise
       call point_to_hydro_grid(mla,variables,p_v, p_rho, p_c, p_T, p_s, &
                                 umac,rho,temperature,scalars)
       if(.false.) write(*,*) "associated=", associated(p_v), associated(p_rho), &
                   associated(p_c), associated(p_T), associated(p_s)
                                 
       call updateHydroAnalysisPrimitive (grid_2D, velocity=p_v, &
           density=p_rho, concentration=p_c, temperature=p_T, scalars=p_s)

       !---------------------------------------------------

       ! Point pointers to the right places, if present, NULL otherwise
       call point_to_hydro_grid(mla,variables_1D,p_v, p_rho, p_c, p_T, p_s, &
                                umac,rho,temperature,scalars)

       call updateHydroAnalysisPrimitive (grid_1D, velocity=p_v, &
           density=p_rho, concentration=p_c, temperature=p_T, scalars=p_s)
                                      
    end if   

    deallocate(variables_1D)

  contains
  
    subroutine average_1D(variables_1D)
    
      real(dp_t), dimension(1:ncells(pdim),nvar) :: variables_1D
    
      integer :: ii, jj, kk
      integer :: lo(3), hi(3) ! It is better to make this dimension independent this way

      variables_1D = 0.0_dp_t
      
      do i=1,nfabs(s_dir)
         sdp => dataptr(s_dir, i)
         lo=1; lo(1:dm) = lwb(get_box(s_dir, i))
         hi=1; hi(1:dm) = upb(get_box(s_dir, i))
         select case(pdim)
         case (1)
            do kk=lo(3),hi(3)
               do jj=lo(2),hi(2)
                  variables_1D = variables_1D + sdp(:,jj,kk,1:nvar)
               end do
            end do
         case (2)
            do kk=lo(3),hi(3)
               do ii=lo(1),hi(1)
                  variables_1D = variables_1D + sdp(ii,:,kk,1:nvar)
               end do
            end do
         case (3)
            do jj=lo(2),hi(2)
               do ii=lo(1),hi(1)
                  variables_1D = variables_1D + sdp(ii,jj,:,1:nvar)
               end do
            end do
         end select
      end do

       ! Divide by number of cells to get average
       ! Note that there is no need to divide again after the MPI reduction
       lo=1; lo(1:dm) = lwb(s_dir%la%lap%pd)
       hi=1; hi(1:dm) = upb(s_dir%la%lap%pd)
       variables_1D = variables_1D / ( (hi(qdim)-lo(qdim)+1)*(hi(qdim2)-lo(qdim2)+1) )
    
    end subroutine 

  end subroutine

  ! This subroutine writes out the mean and variance of all variables using
  ! different averaging techniques.
  !
  ! We work with m_in and s_in (if analyze_conserved=T) or umac and rho (=F).
  !
  ! First, it writes out "vertical" averages, as defined by the
  ! pdim=abs(project_dir) direction.
  ! Thus, in 3D, it writes out a "2D" plotfile (a 3D plotfile with ncell=1
  ! in the pdim direction) called "vstatXXXXXX".
  ! In 2D, it writes out a 1D text file, also called "vstatXXXXX"
  !
  ! Next, it writes out "horizontal" averages.
  ! Thus, in both 2D and 3D, it writes out a 1D text file called "hstatXXXXXX"
  subroutine print_stats(mla,dx,step,time,umac,rho,temperature,scalars)
    type(ml_layout), intent(in   ) :: mla
    real(dp_t)     , intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: step
    real(kind=dp_t), intent(in   ) :: time
    type(multifab) , intent(in), optional :: umac(:,:)
    type(multifab) , intent(in), dimension(:), optional :: rho,temperature,scalars
  
    ! Local
    integer nlevs,dm,pdim,qdim,qdim2,i,n,comp
    integer lo(3),hi(3)
    integer ii,jj,kk

    ! pointers to access s_dir, s_projected, and s_var multifabs
    real(kind=dp_t), pointer :: sdp(:,:,:,:)
    real(kind=dp_t), pointer :: spp(:,:,:,:)
    real(kind=dp_t), pointer :: svp(:,:,:,:)

    real(kind=dp_t), allocatable :: stats_1d(:,:), stats_1d_proc(:,:)

    character(len=20) :: plotfile_name
    character(len=20) :: variable_names(2*nvar)
    type(multifab) :: plotdata(1)
    integer :: rr(0)
    type(box) :: bx

    nlevs = mla%nlevel
    dm    = mla%dim
    pdim  = abs(project_dir)

    ! qdim is orthogonal to pdim
    if (pdim .eq. 1) then
       qdim  = 2
       qdim2 = 3
    else if (pdim .eq. 2) then
       qdim  = 1
       qdim2 = 3
    else if (pdim .eq. 3) then
       qdim  = 1
       qdim2 = 2
    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    ! COMPUTE AND WRITE OUT VERTICAL AVERAGE/VARIANCE
    !!!!!!!!!!!!!!!!!!!!!!!

    ! Re-distribute the full grid into a grid of "tall skinny boxes"
    ! These boxes are not distributed along project_dim so we can do local analysis easily
    ! -------------------------
    call gather_hydro_grid(mla,s_dir,umac,rho,temperature,scalars, variable_names(1:nvar))

    ! Compute s_projected (average) and s_var (variance)
    ! -------------------------
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       spp => dataptr(s_projected, i)
       svp => dataptr(s_var, i)
       lo(1:dm) = lwb(get_box(s_dir, i))
       hi(1:dm) = upb(get_box(s_dir, i))
       ! first put average, <x>, into s_projected
       ! then put variance, < (x - <x> )^2 > = <x^2> - <x>^2, into s_var
       select case (pdim)
       case (1)
          spp(0,:,:,:)=SUM( sdp, DIM=1 ) / dble(hi(1)-lo(1)+1)
          svp(0,:,:,:)=SUM( sdp**2, DIM=1 ) / dble(hi(1)-lo(1)+1) - spp(0,:,:,:)**2
       case (2)
          spp(:,0,:,:)=SUM( sdp, DIM=2 ) / dble(hi(2)-lo(2)+1)
          svp(:,0,:,:)=SUM( sdp**2, DIM=2 ) / dble(hi(2)-lo(2)+1) - spp(:,0,:,:)**2
       case (3)
          spp(:,:,0,:)=SUM( sdp, DIM=3 ) / dble(hi(3)-lo(3)+1)
          svp(:,:,0,:)=SUM( sdp**2, DIM=3 ) / dble(hi(3)-lo(3)+1) - spp(:,:,0,:)**2
       end select
    end do

    ! For 2D simulations, write the vertical average/variance to a text file
    if (dm .eq. 2) then
       
       ! collect the vertical average/variance in a 1D
       ! array so we can write it out in proper order

       bx = s_projected%la%lap%pd
       lo(1:dm) = lwb(bx)
       hi(1:dm) = upb(bx)

       ! components 1:nvar        will hold the average
       ! components nvar+1:2*nvar will hold the variance
       allocate(stats_1d_proc(lo(qdim):hi(qdim),2*nvar))
       allocate(stats_1d     (lo(qdim):hi(qdim),2*nvar))

       stats_1d_proc = 0.d0
       stats_1d      = 0.d0

       ! collect average and variance
       do i=1,nfabs(s_projected)
          spp => dataptr(s_projected, i)
          svp => dataptr(s_var, i)
          lo(1:dm) = lwb(get_box(s_projected, i))
          hi(1:dm) = upb(get_box(s_projected, i))
          select case(pdim)
          case (1)
             stats_1d_proc(lo(2):hi(2),1:nvar)        = spp(0,lo(2):hi(2),1,1:nvar)
             stats_1d_proc(lo(2):hi(2),nvar+1:2*nvar) = svp(0,lo(2):hi(2),1,1:nvar)
          case (2)
             stats_1d_proc(lo(1):hi(1),1:nvar)        = spp(lo(1):hi(1),0,1,1:nvar)
             stats_1d_proc(lo(1):hi(1),nvar+1:2*nvar) = svp(lo(1):hi(1),0,1,1:nvar)
          end select
       end do

       ! sum reduction
       do n=1,2*nvar
          call parallel_reduce(stats_1d(:,n), stats_1d_proc(:,n), MPI_SUM, &
                               proc=parallel_IOProcessorNode())
       end do

       if ( parallel_IOProcessor() ) then

          ! define the name of the statfile that will be written
          write(unit=plotfile_name,fmt='("vstat",i8.8)') step
          write(*,'(2A)') "Saving vSTAT FILEs to file ", trim(plotfile_name)
          write(*,*)

          lo(1:dm) = lwb(bx)
          hi(1:dm) = upb(bx)
          open(1000, file=trim(plotfile_name), status = "unknown", action = "write")

          write(1000,'(A)', advance="no") "# 1=x, "
          do comp=1, nvar
             write(1000,'(I0,3A)', advance="no") comp+1,"=av(",trim(variable_names(comp)), "), "
          end do
          do comp=1, nvar
             write(1000,'(I0,3A)', advance="no") nvar+comp+1,"=var(",trim(variable_names(comp)), "), "
          end do
          write(1000,*) ! New line

          do i=lo(qdim),hi(qdim)
             write(1000,'(1000(g17.9))') prob_lo(qdim) + (i+0.5d0)*dx(1,qdim), &
                  stats_1d(i,:)
          end do
          
          close(1000)
       end if

       deallocate(stats_1d,stats_1d_proc)
   
    ! For 3D simulations, write the vertical average/variance to a plotfile
    else if (dm .eq. 3) then

       ! components 1:nvar        will hold the average
       ! components nvar+1:2*nvar will hold the variance
       call multifab_build(plotdata(1),s_projected%la,2*nvar,0)

       ! copy s_projected and s_var into plotdata
       call multifab_copy_c(plotdata(1),1     ,s_projected,1,nvar)
       call multifab_copy_c(plotdata(1),nvar+1,s_var      ,1,nvar)

       do comp=1, nvar
          variable_names(nvar+comp)="var-"//trim(variable_names(comp))
       end do

       ! define the name of the plotfile that will be written
       write(unit=plotfile_name,fmt='("vstat",i8.8)') step
       if ( parallel_IOProcessor() ) then
          write(*,'(2A)') "Saving vSTAT FILEs to directory ", trim(plotfile_name)
          write(*,*)
       end if

       ! write the plotfile
       call fabio_ml_multifab_write_d(plotdata, rr, plotfile_name, variable_names, &
                                      plotdata(1)%la%lap%pd, prob_lo, prob_hi, &
                                      time, dx(1,:))

       call multifab_destroy(plotdata(1))

    end if

    !!!!!!!!!!!!!!!!!!!!!!!
    ! COMPUTE AND WRITE OUT HORIZONTAL AVERAGE/VARIANCE
    !!!!!!!!!!!!!!!!!!!!!!!

    bx = s_dir%la%lap%pd
    lo=1; lo(1:dm) = lwb(bx)
    hi=1; hi(1:dm) = upb(bx)

    ! will hold average (components 1:nvar) and variance (nvar+1:2*nvar)
    ! as a function of pdim
    allocate(stats_1d_proc(lo(pdim):hi(pdim),2*nvar))
    allocate(stats_1d     (lo(pdim):hi(pdim),2*nvar))

    stats_1d_proc = 0.d0
    stats_1d      = 0.d0

    ! put sum of x into average 
    ! put sum of x^2 and variance
    do i=1,nfabs(s_dir)
       sdp => dataptr(s_dir, i)
       lo=1; lo(1:dm) = lwb(get_box(s_dir, i))
       hi=1; hi(1:dm) = upb(get_box(s_dir, i))
       select case(pdim)
       case (1)
          do kk=lo(3),hi(3)
             do jj=lo(2),hi(2)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(:,jj,kk,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(:,jj,kk,1:nvar)**2
             end do
          end do
       case (2)
          do kk=lo(3),hi(3)
             do ii=lo(1),hi(1)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(ii,:,kk,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(ii,:,kk,1:nvar)**2
             end do
          end do
       case (3)
          do jj=lo(2),hi(2)
             do ii=lo(1),hi(1)
                stats_1d_proc(:,1:nvar) = stats_1d_proc(:,1:nvar) + sdp(ii,jj,:,1:nvar)
                stats_1d_proc(:,nvar+1:2*nvar) = stats_1d_proc(:,nvar+1:2*nvar) &
                     + sdp(ii,jj,:,1:nvar)**2
             end do
          end do
       end select
    end do

    ! sum reduction
    do n=1,2*nvar
       call parallel_reduce(stats_1d(:,n), stats_1d_proc(:,n), MPI_SUM, &
                            proc=parallel_IOProcessorNode())
    end do

    if ( parallel_IOProcessor() ) then

       ! divide by number of cells, so now average contains <x> and variance contains <x^2>
       lo=1; lo(1:dm) = lwb(bx)
       hi=1; hi(1:dm) = upb(bx)
       stats_1d = stats_1d / ( (hi(qdim)-lo(qdim)+1)*(hi(qdim2)-lo(qdim2)+1) )

       ! calculate variance = <x^2> - <x>^2
       stats_1d(:,nvar+1:2*nvar) = stats_1d(:,nvar+1:2*nvar) - stats_1d(:,1:nvar)**2

       ! define the name of the statfile that will be written
       write(unit=plotfile_name,fmt='("hstat",i8.8)') step
       write(*,'(2A)') "Saving hSTAT FILEs to file ", trim(plotfile_name)
       write(*,*)

       lo(1:dm) = lwb(bx)
       hi(1:dm) = upb(bx)
       open(1000, file=trim(plotfile_name), status = "unknown", action = "write")

       write(1000,*) '# time',time

       write(1000,'(A)', advance="no") "# 1=y, "
       do comp=1, nvar
          write(1000,'(I0,3A)', advance="no") comp+1,"=av(",trim(variable_names(comp)), "), "
       end do
       do comp=1, nvar
          write(1000,'(I0,3A)', advance="no") nvar+comp+1,"=var(",trim(variable_names(comp)), "), "
       end do
       write(1000,*) ! New line

       do i=lo(pdim),hi(pdim)
          write(1000,'(1000(g17.9))') prob_lo(pdim) + (i+0.5d0)*dx(1,pdim), stats_1d(i,:)
       end do

       close(1000)

    end if

  end subroutine print_stats

   ! We need to make centered velocities to use FFTs directly
   ! These have a grid of size (nx,ny,nz) and not (nx+1,ny,nz) etc. like the staggered velocities do
   subroutine StaggeredToCentered(mac_cc,mla,umac)
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in) ::    umac(:,:)
    type(multifab), intent(inout) :: mac_cc(mla%nlevel) 

    integer :: i, n, nlevs, dm, pdim

    nlevs = mla%nlevel
    dm = mla%dim
    
    if(center_snapshots) then
       do i=1,dm
          call average_face_to_cc(mla,umac(:,i),1,mac_cc,i,1)
       end do
    else
       ! Pretend that velocities are cell-centered instead of staggered
       ! This is not right but for periodic should be OK as it preserves the spectrum of fluctuations
       do i=1,dm
          call shift_face_to_cc(mla,umac(:,i),1,mac_cc,i,1)
       end do       
    end if
  
  end subroutine  

end module
