module restart_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use checkpoint_module
  use define_bc_module
  use probin_common_module, only: dim_in, advection_type, restart, advection_type, algorithm_type
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: initialize_from_restart

contains

  subroutine initialize_from_restart(mla,time,dt,rho,rhotot,pres,diff_mass_fluxdiv, &
                                     stoch_mass_fluxdiv,umac,pmask)
 
     type(ml_layout),intent(out)   :: mla
     real(dp_t)    , intent(  out) :: time,dt
     type(multifab), intent(inout) :: rho(:)
     type(multifab), intent(inout) :: rhotot(:)
     type(multifab), intent(inout) :: pres(:)
     type(multifab), intent(inout) :: diff_mass_fluxdiv(:)
     type(multifab), intent(inout) :: stoch_mass_fluxdiv(:)
     type(multifab), intent(inout) :: umac(:,:)
     logical       , intent(in   ) :: pmask(:)

     type(ml_boxarray)         :: mba
     type(multifab), pointer   :: chkdata(:)
     type(multifab), pointer   :: chkdata_edgex(:)
     type(multifab), pointer   :: chkdata_edgey(:)
     type(multifab), pointer   :: chkdata_edgez(:)
     type(layout)              :: la

     integer :: n,nlevs,i,dm,ng_s

     type(bl_prof_timer), save :: bpt

     call build(bpt,"initialize_from_restart")

     dm = dim_in

     if (advection_type .eq. 0) then
        ng_s = 2 ! centered advection
     else if (advection_type .le. 3) then
        ng_s = 3 ! bilinear bds or unlimited quadratic bds
     else if (advection_type .eq. 4) then
        ng_s = 4 ! limited quadratic bds
     end if

     call fill_restart_data(mba,chkdata,chkdata_edgex,chkdata_edgey,chkdata_edgez, &
                            time,dt)

     call ml_layout_build(mla,mba,pmask)

     nlevs = mba%nlevel

     do n = 1,nlevs
        call multifab_build(rho(n)   , mla%la(n), nspecies, ng_s)
        call multifab_build(rhotot(n), mla%la(n),        1, ng_s)
        call multifab_build(pres(n)  , mla%la(n),        1, 1)
        call multifab_build(diff_mass_fluxdiv(n), mla%la(n),nspecies,0) 
        call multifab_build(stoch_mass_fluxdiv(n),mla%la(n),nspecies,0) 
        do i=1,dm
           call multifab_build_edge(umac(n,i), mla%la(n), 1, 1, i)
        end do
     end do
     do n = 1,nlevs
        call multifab_copy_c(rho(n)   ,1,chkdata(n)      ,1         ,nspecies)
        call multifab_copy_c(rhotot(n),1,chkdata(n)      ,nspecies+1,1)
        call multifab_copy_c(pres(n)  ,1,chkdata(n)      ,nspecies+2,1)
        if (algorithm_type .eq. 0) then
           call multifab_copy_c( diff_mass_fluxdiv(n),1,chkdata(n),  nspecies+3,nspecies)
           call multifab_copy_c(stoch_mass_fluxdiv(n),1,chkdata(n),2*nspecies+3,nspecies)
        end if
        call multifab_copy_c(umac(n,1),1,chkdata_edgex(n),1         ,1)
        call multifab_copy_c(umac(n,2),1,chkdata_edgey(n),1         ,1)
        if (dm .eq. 3) then
           call multifab_copy_c(umac(n,3),1,chkdata_edgez(n),1,1)
        end if
        !
        ! The layout for chkdata is built standalone, level
        ! by level, and need to be destroy()d as such as well.
        !
        la = get_layout(chkdata(n))
        call multifab_destroy(chkdata(n))
        call destroy(la)
        la = get_layout(chkdata_edgex(n))
        call multifab_destroy(chkdata_edgex(n))
        call destroy(la)
        la = get_layout(chkdata_edgey(n))
        call multifab_destroy(chkdata_edgey(n))
        call destroy(la)
        if (dm .eq. 3) then
           la = get_layout(chkdata_edgez(n))
           call multifab_destroy(chkdata_edgez(n))
           call destroy(la)
        end if
     end do
     deallocate(chkdata)
     deallocate(chkdata_edgex)
     deallocate(chkdata_edgey)
     if (dm .eq. 3) then
        deallocate(chkdata_edgez)
     end if
     call destroy(mba)

     call destroy(bpt)

  end subroutine initialize_from_restart

  subroutine fill_restart_data(mba,chkdata, &
                               chkdata_edgex,chkdata_edgey,chkdata_edgez, &
                               time,dt)


    real(dp_t)       , intent(  out) :: time,dt
    type(ml_boxarray), intent(  out) :: mba

    type(multifab)   , pointer        :: chkdata(:)
    type(multifab)   , pointer        :: chkdata_edgex(:)
    type(multifab)   , pointer        :: chkdata_edgey(:)
    type(multifab)   , pointer        :: chkdata_edgez(:)
    character(len=11)                 :: sd_name
    integer                           :: n,nlevs,dm
    integer                           :: rrs(10)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_restart_data")

    dm = dim_in

    write(unit=sd_name,fmt='("chk",i8.8)') restart

    if ( parallel_IOProcessor() ) then
       print *,'Reading ',sd_name,' to get state data for restart'
    end if

    call checkpoint_read(chkdata, chkdata_edgex, chkdata_edgey, chkdata_edgez, &
                         sd_name, rrs, time, dt, nlevs)

    call build(mba,nlevs,dm)
    mba%pd(1) =  bbox(get_boxarray(chkdata(1)))
    do n = 2,nlevs
      mba%pd(n) = refine(mba%pd(n-1),2)
      mba%rr(n-1,:) = rrs(n-1)
    end do
    do n = 1,nlevs
      call boxarray_build_copy(mba%bas(n), get_boxarray(chkdata(n))) 
    end do

    call destroy(bpt)

  end subroutine fill_restart_data

end module restart_module
