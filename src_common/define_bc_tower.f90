module define_bc_module

  use bl_types
  use ml_layout_module
  use bc_module
  use inhomogeneous_bc_val_module

  implicit none

  type bc_level

     integer, pointer :: phys_bc_level_array(:,:,:) => Null()
     integer, pointer ::  adv_bc_level_array(:,:,:,:) => Null()
     integer, pointer ::  ell_bc_level_array(:,:,:,:) => Null()

  end type bc_level

  type bc_tower

     integer :: max_level_built = 0
     type(bc_level), pointer :: bc_tower_array(:) => Null()
     integer       , pointer :: domain_bc(:,:) => Null()

  end type bc_tower

  private

  public :: bc_level, bc_tower, &
            initialize_bc, bc_tower_init, bc_tower_level_build, bc_tower_destroy

contains

  subroutine initialize_bc(the_bc_tower,num_levs,dm,pmask,num_scal_bc_in,num_tran_bc_in)

     use bc_module

     use probin_common_module, only : bc_lo, bc_hi

     type(bc_tower), intent(  out) :: the_bc_tower
     integer       , intent(in   ) :: num_levs,dm,num_scal_bc_in,num_tran_bc_in
     logical       , intent(in   ) :: pmask(:)

     integer :: domain_phys_bc(dm,2)

     num_scal_bc  = num_scal_bc_in
     num_tran_bc  = num_tran_bc_in

     vel_bc_comp  = 1
     pres_bc_comp = vel_bc_comp+dm
     scal_bc_comp = pres_bc_comp+1
     tran_bc_comp = scal_bc_comp+num_scal_bc

     ! Define the physical boundary conditions on the domain
     ! Put the bc values from the inputs file into domain_phys_bc
     domain_phys_bc(1,1) = bc_lo(1)
     domain_phys_bc(1,2) = bc_hi(1)
     if (pmask(1)) then
        domain_phys_bc(1,:) = BC_PER
        if (bc_lo(1) .ne. -1 .or. bc_hi(1) .ne. -1) &
             call bl_error('MUST HAVE BCX = -1 if PMASK = T')
     end if
     if (dm > 1) then
        domain_phys_bc(2,1) = bc_lo(2)
        domain_phys_bc(2,2) = bc_hi(2)
        if (pmask(2)) then
           domain_phys_bc(2,:) = BC_PER
           if (bc_lo(2) .ne. -1 .or. bc_hi(2) .ne. -1) &
                call bl_error('MUST HAVE BCY = -1 if PMASK = T') 
        end if
     end if
     if (dm > 2) then
        domain_phys_bc(3,1) = bc_lo(3)
        domain_phys_bc(3,2) = bc_hi(3)
        if (pmask(3)) then
           domain_phys_bc(3,:) = BC_PER
           if (bc_lo(3) .ne. -1 .or. bc_hi(3) .ne. -1) &
                call bl_error('MUST HAVE BCZ = -1 if PMASK = T')
        end if
     end if

     ! Initialize the_bc_tower object.
     call bc_tower_init(the_bc_tower,num_levs,dm,domain_phys_bc)

  end subroutine initialize_bc

  subroutine bc_tower_init(bct,num_levs,dm,phys_bc_in)

    type(bc_tower ), intent(  out) :: bct
    integer        , intent(in   ) :: num_levs
    integer        , intent(in   ) :: dm
    integer        , intent(in   ) :: phys_bc_in(:,:)

    allocate(bct%bc_tower_array(num_levs))
    allocate(bct%domain_bc(dm,2))

    bct%domain_bc(:,:) = phys_bc_in(:,:)

  end subroutine bc_tower_init

  subroutine bc_tower_level_build(bct,n,la)

    type(bc_tower ), intent(inout) :: bct
    integer        , intent(in   ) :: n
    type(layout)   , intent(in   ) :: la

    integer :: ngrids,dm
    integer :: default_value

    if (associated(bct%bc_tower_array(n)%phys_bc_level_array)) then
      deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
      deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
      deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
      bct%bc_tower_array(n)%phys_bc_level_array => NULL()
      bct%bc_tower_array(n)%adv_bc_level_array => NULL()
      bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end if

    ngrids = layout_nlocal(la)
    dm = layout_dim(la)

    allocate(bct%bc_tower_array(n)%phys_bc_level_array(0:ngrids,dm,2))
    default_value = INTERIOR
    call phys_bc_level_build(bct%bc_tower_array(n)%phys_bc_level_array,la, &
                             bct%domain_bc,default_value)

    ! Here we allocate (dm+1)=pres_bc_comp components for x_u and x_p
    !                  num_scal_bc components for scalars
    !                  num_tran_bc components for transport coefficients
    allocate(bct%bc_tower_array(n)%adv_bc_level_array(0:ngrids,dm,2,pres_bc_comp+num_scal_bc+num_tran_bc))
    default_value = INTERIOR
    call adv_bc_level_build(bct%bc_tower_array(n)%adv_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

    ! Here we allocate (dm+1)=pres_bc_comp components for x_u and x_p
    !                  num_scal_bc components for scalars
    !                  num_tran_bc components for transport coefficients
    allocate(bct%bc_tower_array(n)%ell_bc_level_array(0:ngrids,dm,2,pres_bc_comp+num_scal_bc+num_tran_bc))
    default_value = BC_INT
    call ell_bc_level_build(bct%bc_tower_array(n)%ell_bc_level_array, &
                            bct%bc_tower_array(n)%phys_bc_level_array,default_value)

     bct%max_level_built = n

  end subroutine bc_tower_level_build

  subroutine bc_tower_destroy(bct)

    type(bc_tower), intent(inout) :: bct

    integer :: n

    do n = 1,bct%max_level_built
       deallocate(bct%bc_tower_array(n)%phys_bc_level_array)
       deallocate(bct%bc_tower_array(n)%adv_bc_level_array)
       deallocate(bct%bc_tower_array(n)%ell_bc_level_array)
       bct%bc_tower_array(n)%phys_bc_level_array => NULL()
       bct%bc_tower_array(n)%adv_bc_level_array => NULL()
       bct%bc_tower_array(n)%ell_bc_level_array => NULL()
    end do
    deallocate(bct%bc_tower_array)

    deallocate(bct%domain_bc)

  end subroutine bc_tower_destroy

  subroutine phys_bc_level_build(phys_bc_level,la_level,domain_bc,default_value)

    integer     , intent(inout) :: phys_bc_level(0:,:,:)
    integer     , intent(in   ) :: domain_bc(:,:)
    type(layout), intent(in   ) :: la_level
    integer     , intent(in   ) :: default_value
    type(box) :: bx,pd
    integer :: d,i

    pd = layout_get_pd(la_level) 

    phys_bc_level = default_value

    i = 0
    do d = 1,layout_dim(la_level)
       phys_bc_level(i,d,1) = domain_bc(d,1)
       phys_bc_level(i,d,2) = domain_bc(d,2)
    end do

    do i = 1,layout_nlocal(la_level)
       bx = layout_get_box(la_level,global_index(la_level,i))
       do d = 1,layout_dim(la_level)
          if (lwb(bx,d) == lwb(pd,d)) phys_bc_level(i,d,1) = domain_bc(d,1)
          if (upb(bx,d) == upb(pd,d)) phys_bc_level(i,d,2) = domain_bc(d,2)
       end do
    end do

  end subroutine phys_bc_level_build

  subroutine adv_bc_level_build(adv_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  adv_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm,igrid,d,lohi

    ! these boundary conditions are written for a stokes system where we 
    ! solve for the INCREMENT to velocity and pressure

    adv_bc_level = default_value

    dm = size(adv_bc_level,dim=2)

    do igrid = 0, size(adv_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2

       if (phys_bc_level(igrid,d,lohi) == PERIODIC .or. &
           phys_bc_level(igrid,d,lohi) == INTERIOR ) then

          ! retain the default value of INTERIOR

       else if ( (phys_bc_level(igrid,d,lohi) >= NO_SLIP_START) .and.  (phys_bc_level(igrid,d,lohi) <= NO_SLIP_END) ) then

          ! for normal velocity we impose a Dirichlet velocity condition
          ! for transverse velocity we imposie a Dirichlet velocity condition
          ! for pressure we use extrapolation
          adv_bc_level(igrid,d,lohi,1:dm)                    = DIR_VEL  ! normal and transverse velocity
          adv_bc_level(igrid,d,lohi,pres_bc_comp)                    = FOEXTRAP ! pressure

       else if ( (phys_bc_level(igrid,d,lohi) >= SLIP_START) .and.  (phys_bc_level(igrid,d,lohi) <= SLIP_END) ) then

          ! for normal velocity we impose a Dirichlet velocity condition
          ! for transverse velocity we imposie a Dirichlet stress condition
          ! for pressure we use extrapolation
          adv_bc_level(igrid,d,lohi,vel_bc_comp:vel_bc_comp+dm-1) = DIR_TRACT ! transverse velocity
          adv_bc_level(igrid,d,lohi,vel_bc_comp+d-1)              = DIR_VEL   ! normal velocity
          adv_bc_level(igrid,d,lohi,pres_bc_comp)                 = FOEXTRAP  ! pressure

       else

          print*,'adv_bc_level_build',igrid,d,lohi,phys_bc_level(igrid,d,lohi)
          call bl_error('BC TYPE NOT SUPPORTED 1')

       end if
       
       ! The scalars can be handled in many different ways, so defer to application-specific code:
       call scalar_bc(phys_bc_level(igrid,d,lohi),adv_bc_level(igrid,d,lohi, &
                      scal_bc_comp:scal_bc_comp+num_scal_bc-1))
       
       ! The transport coefficients can be handled in many different ways, so defer to application-specific code:
       call transport_bc(phys_bc_level(igrid,d,lohi),adv_bc_level(igrid,d,lohi, &
                         tran_bc_comp:tran_bc_comp+num_tran_bc-1))

       if (any(adv_bc_level(igrid,d,lohi,scal_bc_comp:scal_bc_comp+num_scal_bc-1) .eq. -999)) then

          print*,'adv_bc_level_build',igrid,d,lohi,phys_bc_level(igrid,d,lohi)
          call bl_error('BC TYPE NOT SUPPORTED 2')

       end if

    end do
    end do
    end do

  end subroutine adv_bc_level_build

  subroutine ell_bc_level_build(ell_bc_level,phys_bc_level,default_value)

    integer  , intent(inout) ::  ell_bc_level(0:,:,:,:)
    integer  , intent(in   ) :: phys_bc_level(0:,:,:)
    integer  , intent(in   ) :: default_value

    integer :: dm
    integer :: igrid,d,lohi

    ell_bc_level = default_value
 
    dm = size(ell_bc_level,dim=2)

    do igrid = 0, size(ell_bc_level,dim=1)-1
    do d = 1, dm
    do lohi = 1, 2

       if (phys_bc_level(igrid,d,lohi) == INTERIOR) then

          ! retain the default value of INTERIOR

       else if (phys_bc_level(igrid,d,lohi) == PERIODIC) then

          ! pressure is periodic
          ell_bc_level(igrid,d,lohi,pres_bc_comp) = BC_PER

       else

          ! pressure and temperature are homogeneous neumann
          ell_bc_level(igrid,d,lohi,pres_bc_comp) = BC_NEU

       end if

    end do
    end do
    end do

  end subroutine ell_bc_level_build

end module define_bc_module
