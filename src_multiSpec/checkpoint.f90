module checkpoint_module

  use parallel
  use bl_types
  use multifab_module
  use ml_layout_module
  use bl_IO_module
  use fab_module
  use fabio_module, only: fabio_mkdir, fabio_ml_multifab_write_d
  use probin_common_module, only: dim_in, algorithm_type
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: checkpoint_write, checkpoint_read

contains

  subroutine checkpoint_write(mla,rho,rhotot,pres,diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                              umac,time,dt,istep_to_write)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)                ! cell-centered partial densities
    type(multifab) , intent(in   ) :: rhotot(:)             ! cell-centered total density
    type(multifab) , intent(in   ) :: pres(:)               ! cell-centered pressure
    type(multifab) , intent(in   ) :: diff_mass_fluxdiv(:)  ! diffusive mass fluxes
    type(multifab) , intent(in   ) :: stoch_mass_fluxdiv(:) ! stochastic mass fluxes
    type(multifab) , intent(in   ) :: umac(:,:)             ! edge-based velocities
    integer        , intent(in   ) :: istep_to_write
    real(kind=dp_t), intent(in   ) :: time,dt

    type(multifab), pointer :: chkdata(:)
    type(multifab), pointer :: chkdata_edge(:,:)

    integer :: n,nlevs,dm,i

    character(len=11) :: sd_name

    type(bl_prof_timer), save :: bpt

    call build(bpt, "checkpoint_write")

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(chkdata(nlevs))
    allocate(chkdata_edge(nlevs,dm))
    do n = 1,nlevs
       if (algorithm_type .eq. 0) then
          ! nspecies densities + 1 total density + 1 pressure + diff/stoch_mass_fluxdiv
          call multifab_build(chkdata(n), mla%la(n), 3*nspecies+2, 0)
       else
          ! nspecies densities + 1 total density + 1 pressure
          call multifab_build(chkdata(n), mla%la(n), nspecies+2, 0)
       end if
       ! copy partial densities, total density, and pressure
       call multifab_copy_c(chkdata(n), 1         , rho(n)   , 1, nspecies)
       call multifab_copy_c(chkdata(n), nspecies+1, rhotot(n), 1, 1)
       call multifab_copy_c(chkdata(n), nspecies+2, pres(n)  , 1, 1)
       if (algorithm_type .eq. 0) then
          call multifab_copy_c(chkdata(n), nspecies+3, diff_mass_fluxdiv(n), 1, nspecies)
          call multifab_copy_c(chkdata(n), 2*nspecies+3, stoch_mass_fluxdiv(n), 1, nspecies)
       end if
       do i=1,dm
          ! 1 velocity component and 1 normal bc component for each face
          call multifab_build_edge(chkdata_edge(n,i), mla%la(n), 1, 0, i)
          ! copy velocities
          call multifab_copy_c(chkdata_edge(n,i), 1, umac(n,i), 1, 1)
       end do
    end do
    write(unit=sd_name,fmt='("chk",i8.8)') istep_to_write

    call checkpoint_write_doit(nlevs, sd_name, chkdata, chkdata_edge, mla%mba%rr, time, dt)

    do n = 1,nlevs
       call multifab_destroy(chkdata(n))
       do i=1,dm
          call multifab_destroy(chkdata_edge(n,i))
       end do
    end do
    deallocate(chkdata)
    deallocate(chkdata_edge)

    call destroy(bpt)

  contains

    subroutine checkpoint_write_doit(nlevs_in, dirname, mfs, mfs_edge, rrs, time_in, dt_in)
      
      integer         , intent(in) :: nlevs_in
      type(multifab)  , intent(in) :: mfs(:)
      type(multifab)  , intent(in) :: mfs_edge(:,:)
      integer         , intent(in) :: rrs(:,:)
      character(len=*), intent(in) :: dirname
      real(kind=dp_t) , intent(in) :: time_in, dt_in

      integer :: n
      character(len=128) :: header, sd_name
      integer :: un, dm

      integer         :: nlevs
      real(kind=dp_t) :: time, dt

      type(bl_prof_timer), save :: bpt

      namelist /chkpoint/ time
      namelist /chkpoint/ dt
      namelist /chkpoint/ nlevs

      call build(bpt,"checkpoint_write_doit")

      dm = dim_in
      if ( parallel_IOProcessor() ) call fabio_mkdir(dirname)

      call parallel_barrier() ! All CPUs have to wait till the directory is built.

      write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs, rrs(:,1), sd_name)

      if (parallel_IOProcessor()) then
         print *,'Writing cc state to checkpoint file ',trim(sd_name)
         print *,' '
      end if

      write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_edge(:,1), rrs(:,1), sd_name)

      write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
      call fabio_ml_multifab_write_d(mfs_edge(:,2), rrs(:,1), sd_name)

      if (dm .eq. 3) then
         write(unit=sd_name, fmt='(a,"/State_edgez")') trim(dirname)
         call fabio_ml_multifab_write_d(mfs_edge(:,3), rrs(:,1), sd_name)
      end if

      time  = time_in
      dt = dt_in
      nlevs = nlevs_in

      if (parallel_IOProcessor()) then
         header = "Header"
         un = unit_new()
         open(unit=un, &
              file = trim(dirname) // "/" // trim(header), &
              form = "formatted", access = "sequential", &
              status = "replace", action = "write")
         write(unit=un, nml = chkpoint)
         do n = 1,nlevs-1
            write(unit=un,fmt=*) rrs(n,1)
         end do
         close(un)
      end if

      call destroy(bpt)

    end subroutine checkpoint_write_doit

  end subroutine checkpoint_write

  subroutine checkpoint_read(mfs, mfs_edgex, mfs_edgey, mfs_edgez, &
                             dirname, rrs_out, time_out, dt_out, nlevs_out)

    use bl_IO_module
    use fab_module
    use fabio_module, only: fabio_ml_multifab_read_d
    use parallel

    type(multifab)  ,                pointer :: mfs(:)
    type(multifab)  ,                pointer :: mfs_edgex(:)
    type(multifab)  ,                pointer :: mfs_edgey(:)
    type(multifab)  ,                pointer :: mfs_edgez(:)
    character(len=*), intent(in   )          :: dirname
    integer         , intent(  out)          :: nlevs_out
    real(kind=dp_t) , intent(  out)          :: time_out, dt_out
    integer         , intent(  out)          :: rrs_out(:)

    character(len=128) :: header, sd_name

    integer :: n, un, nlevs, dm
    integer, pointer :: rrs(:)

    real(kind=dp_t) :: time, dt

    type(bl_prof_timer), save :: bpt

    namelist /chkpoint/ nlevs
    namelist /chkpoint/ time
    namelist /chkpoint/ dt

    call build(bpt,"checkpoint_read")

!   First read the header information
    header = "Header"
    un = unit_new()
    open(unit=un, &
         file = trim(dirname) // "/" // trim(header), &
         status = "old", &
         action = "read")
    read(unit=un, nml = chkpoint)
    allocate(rrs(nlevs-1))
    do n = 1,nlevs-1
       read(unit=un,fmt=*) rrs(n)
    end do
    close(un)

     time_out = time
       dt_out = dt
    nlevs_out = nlevs

!   Read the state data into a multilevel multifab.
    write(unit=sd_name, fmt='(a,"/State")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs, sd_name)

    dm = get_dim(mfs(1))

    write(unit=sd_name, fmt='(a,"/State_edgex")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgex, sd_name)

    write(unit=sd_name, fmt='(a,"/State_edgey")') trim(dirname)
    call fabio_ml_multifab_read_d(mfs_edgey, sd_name)

    if (dm .eq. 3) then
       write(unit=sd_name, fmt='(a,"/State_edgez")') trim(dirname)
       call fabio_ml_multifab_read_d(mfs_edgez, sd_name)
    end if

    rrs_out(1:nlevs-1) = rrs(1:nlevs-1)

    deallocate(rrs)

    call destroy(bpt)

  end subroutine checkpoint_read

end module checkpoint_module
