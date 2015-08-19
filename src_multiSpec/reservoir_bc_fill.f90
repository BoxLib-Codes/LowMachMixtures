module reservoir_bc_fill_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: rhobar

  implicit none

  private

  public :: reservoir_bc_fill

contains

  subroutine reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: flux_total(:,:) ! should contain sum of diffusive and stochastic mass fluxes
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: i,dm,n,nlevs,ng_f,ng_b
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fpx(:,:,:,:)
    real(kind=dp_t), pointer :: fpy(:,:,:,:)
    real(kind=dp_t), pointer :: fpz(:,:,:,:)
    real(kind=dp_t), pointer :: vpx(:,:,:,:)
    real(kind=dp_t), pointer :: vpy(:,:,:,:)
    real(kind=dp_t), pointer :: vpz(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "reservoir_bc_fill")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_f = flux_total(1,1)%ng
    ng_b = vel_bc_n(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(flux_total(1,1))
          fpx => dataptr(flux_total(n,1), i)
          fpy => dataptr(flux_total(n,2), i)
          vpx => dataptr(vel_bc_n(n,1), i)
          vpy => dataptr(vel_bc_n(n,2), i)
         lo =  lwb(get_box(flux_total(n,1), i))
         hi =  upb(get_box(flux_total(n,1), i))
          select case (dm)
          case (2)
             call reservoir_bc_fill_2d(fpx(:,:,1,:),fpy(:,:,1,:),ng_f, &
                                       vpx(:,:,1,1),vpy(:,:,1,1),ng_b, &
                                       lo,hi,the_bc_level(n)%phys_bc_level_array(i,:,:))
          case (3)
             fpz => dataptr(flux_total(n,3), i)
             vpz => dataptr(vel_bc_n(n,3), i)
             call reservoir_bc_fill_3d(fpx(:,:,:,:),fpy(:,:,:,:),fpz(:,:,:,:),ng_f, &
                                       vpx(:,:,:,1),vpy(:,:,:,1),vpz(:,:,:,1),ng_b, &
                                       lo,hi,the_bc_level(n)%phys_bc_level_array(i,:,:))
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine reservoir_bc_fill

  subroutine reservoir_bc_fill_2d(fluxx,fluxy,ng_f,vel_bc_nx,vel_bc_ny,ng_b,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_b
    real(kind=dp_t), intent(in   ) ::     fluxx(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(in   ) ::     fluxy(lo(1)-ng_f:,lo(2)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:)
    real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j,n
    real(kind=dp_t) :: sum

    ! update umac bc based on diffusive flux at boundary
    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxx(lo(1),j,n)/rhobar(n)
          end do
          vel_bc_nx(lo(1),j) = sum
       end do
    end if
    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxx(hi(1)+1,j,n)/rhobar(n)
          end do
          vel_bc_nx(hi(1)+1,j) = sum
       end do
    end if

    ! update vmac bc based on diffusive flux at boundary
    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxy(i,lo(2),n)/rhobar(n)
          end do
          vel_bc_ny(i,lo(2)) = sum
       end do
    end if
    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxy(i,hi(2)+1,n)/rhobar(n)
          end do
          vel_bc_ny(i,hi(2)+1) = sum
       end do
    end if

  end subroutine reservoir_bc_fill_2d

  subroutine reservoir_bc_fill_3d(fluxx,fluxy,fluxz,ng_f,vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_b,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_b
    real(kind=dp_t), intent(in   ) ::     fluxx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) ::     fluxy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(in   ) ::     fluxz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)
    real(kind=dp_t), intent(inout) :: vel_bc_nx(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(inout) :: vel_bc_ny(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    real(kind=dp_t), intent(inout) :: vel_bc_nz(lo(1)-ng_b:,lo(2)-ng_b:,lo(3)-ng_b:)
    integer        , intent(in   ) :: bc(:,:)

    ! local
    integer :: i,j,k,n
    real(kind=dp_t) :: sum

    ! update umac bc based on diffusive flux at boundary
    if (bc(1,1) .eq. NO_SLIP_RESERVOIR .or. bc(1,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxx(lo(1),j,k,n)/rhobar(n)
          end do
          vel_bc_nx(lo(1),j,k) = sum
       end do
       end do
    end if
    if (bc(1,2) .eq. NO_SLIP_RESERVOIR .or. bc(1,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxx(hi(1)+1,j,k,n)/rhobar(n)
          end do
          vel_bc_nx(hi(1)+1,j,k) = sum
       end do
       end do
    end if

    ! update vmac bc based on diffusive flux at boundary
    if (bc(2,1) .eq. NO_SLIP_RESERVOIR .or. bc(2,1) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxy(i,lo(2),k,n)/rhobar(n)
          end do
          vel_bc_ny(i,lo(2),k) = sum
       end do
       end do
    end if
    if (bc(2,2) .eq. NO_SLIP_RESERVOIR .or. bc(2,2) .eq. SLIP_RESERVOIR) then
       do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxy(i,hi(2)+1,k,n)/rhobar(n)
          end do
          vel_bc_ny(i,hi(2)+1,k) = sum
       end do
       end do
    end if

    ! update wmac bc based on diffusive flux at boundary
    if (bc(3,1) .eq. NO_SLIP_RESERVOIR .or. bc(3,1) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxz(i,j,lo(3),n)/rhobar(n)
          end do
          vel_bc_nz(i,j,lo(3)) = sum
       end do
       end do
    end if
    if (bc(3,2) .eq. NO_SLIP_RESERVOIR .or. bc(3,2) .eq. SLIP_RESERVOIR) then
       do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          sum = 0.d0
          do n=1,nspecies
             sum = sum + fluxz(i,j,hi(3)+1,n)/rhobar(n)
          end do
          vel_bc_nz(i,j,hi(3)+1) = sum
       end do
       end do
    end if

  end subroutine reservoir_bc_fill_3d

end module reservoir_bc_fill_module
