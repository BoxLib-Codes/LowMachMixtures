module fill_rho_ghost_cells_module

  use multifab_module
  use define_bc_module
  use bc_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: fill_rho_ghost_cells

contains

  subroutine fill_rho_ghost_cells(conc,rhotot,the_bc_level)

    type(multifab) , intent(in   ) :: conc
    type(multifab) , intent(inout) :: rhotot
    type(bc_level) , intent(in   ) :: the_bc_level

    ! local
    integer :: i,dm,ng_c,ng_r
    integer :: lo(get_dim(conc)),hi(get_dim(conc))
    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: rp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_rho_ghost_cells")

    dm = get_dim(conc)
    ng_c = conc%ng
    ng_r = rhotot%ng

    call multifab_fill_boundary(rhotot)

    do i=1,nfabs(conc)
       pp => dataptr(conc,i)
       rp => dataptr(rhotot,i)
       lo = lwb(get_box(conc,i))
       hi = upb(get_box(conc,i))
       select case (dm)
       case (2)
          call fill_rho_ghost_cells_2d(pp(:,:,1,:),ng_c,rp(:,:,1,1),ng_r,lo,hi, &
                                       the_bc_level%adv_bc_level_array(i,:,:,scal_bc_comp))
       case (3)
          call fill_rho_ghost_cells_3d(pp(:,:,:,:),ng_c,rp(:,:,:,1),ng_r,lo,hi, &
                                       the_bc_level%adv_bc_level_array(i,:,:,scal_bc_comp))
       end select
    end do

    call destroy(bpt)

  end subroutine fill_rho_ghost_cells

  subroutine fill_rho_ghost_cells_2d(conc,ng_c,rhotot,ng_r,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_r
    real(kind=dp_t), intent(in   ) ::   conc(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t), intent(inout) :: rhotot(lo(1)-ng_r:,lo(2)-ng_r:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j,n
    real(kind=dp_t) :: rhoinv

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_r,hi(2)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(lo(1)-1,j,n)/rhobar(n)
          end do
          rhotot(lo(1)-ng_r:lo(1)-1,j) = 1.d0/rhoinv
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       do j=lo(2)-ng_r,hi(2)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(hi(1)+1,j,n)/rhobar(n)
          end do
          rhotot(hi(1)+1:hi(1)+ng_r,j) = 1.d0/rhoinv
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,lo(2)-1,n)/rhobar(n)
          end do
          rhotot(i,lo(2)-ng_r:lo(2)-1) = 1.d0/rhoinv
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(2,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,hi(2)+1,n)/rhobar(n)
          end do
          rhotot(i,hi(2)+1:hi(2)+ng_r) = 1.d0/rhoinv
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_2d: bc(2,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine fill_rho_ghost_cells_2d

  subroutine fill_rho_ghost_cells_3d(conc,ng_c,rhotot,ng_r,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_r
    real(kind=dp_t), intent(in   ) ::   conc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t), intent(inout) :: rhotot(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j,k,n
    real(kind=dp_t) :: rhoinv


    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_r,hi(3)+ng_r
       do j=lo(2)-ng_r,hi(2)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(lo(1)-1,j,k,n)/rhobar(n)
          end do
          rhotot(lo(1)-ng_r:lo(1)-1,j,k) = 1.d0/rhoinv
       end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
       do k=lo(3)-ng_r,hi(3)+ng_r
       do j=lo(2)-ng_r,hi(2)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(hi(1)+1,j,k,n)/rhobar(n)
          end do
          rhotot(hi(1)+1:hi(1)+ng_r,j,k) = 1.d0/rhoinv
       end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(1,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
       do k=lo(3)-ng_r,hi(3)+ng_r
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,lo(2)-1,k,n)/rhobar(n)
          end do
          rhotot(i,lo(2)-ng_r:lo(2)-1,k) = 1.d0/rhoinv
       end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(2,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
       do k=lo(3)-ng_r,hi(3)+ng_r
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,hi(2)+1,k,n)/rhobar(n)
          end do
          rhotot(i,hi(2)+1:hi(2)+ng_r,k) = 1.d0/rhoinv
       end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(2,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. EXT_DIR) then
       do j=lo(2)-ng_r,hi(2)+ng_r
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,j,lo(3)-1,n)/rhobar(n)
          end do
          rhotot(i,j,lo(3)-ng_r:lo(3)-1) = 1.d0/rhoinv
       end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(3,1) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

    if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP .or. bc(3,2) .eq. EXT_DIR) then
       do j=lo(2)-ng_r,hi(2)+ng_r
       do i=lo(1)-ng_r,hi(1)+ng_r
          rhoinv = 0.d0
          do n=1,nspecies
             rhoinv = rhoinv + conc(i,j,hi(3)+1,n)/rhobar(n)
          end do
          rhotot(i,j,hi(3)+1:hi(3)+ng_r) = 1.d0/rhoinv
       end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'fill_rho_ghost_cells_3d: bc(3,2) =',bc(1,1)
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine fill_rho_ghost_cells_3d

end module fill_rho_ghost_cells_module
