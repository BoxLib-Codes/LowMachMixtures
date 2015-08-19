module multifab_coefbc_module

  use multifab_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: multifab_coefbc

contains

  subroutine multifab_coefbc(s,start_scomp,ncomp,the_bc_level)

    ! this fills ghost cells for transport coefficients at walls by
    ! averaging ghost and interior into ghost

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,ncomp
    type(bc_level) , intent(in   ) :: the_bc_level
   
    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng,dm,comp
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"multifab_coefbc")

    ng = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       do comp=start_scomp,start_scomp+ncomp-1
          select case (dm)
          case (2)
             call coefbc_2d(sp(:,:,1,comp), lo, hi, ng, &
                            the_bc_level%phys_bc_level_array(i,:,:))
          case (3)
             call coefbc_3d(sp(:,:,:,comp), lo, hi, ng, &
                            the_bc_level%phys_bc_level_array(i,:,:))
          end select
       end do
    end do
 
    call destroy(bpt)

  end subroutine multifab_coefbc

  subroutine coefbc_2d(s,lo,hi,ng,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-ng:lo(1)-1,j) = 0.5d0*(s(lo(1),j)+s(lo(1)-1,j))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1:hi(1)+ng,j) = 0.5d0*(s(hi(1),j)+s(hi(1)+1,j))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-ng:lo(2)-1) = 0.5d0*(s(i,lo(2))+s(i,lo(2)-1))
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1:hi(2)+ng) = 0.5d0*(s(i,hi(2))+s(i,hi(2)+1))
       end do
    end if

  end subroutine coefbc_2d

  subroutine coefbc_3d(s,lo,hi,ng,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)

    ! Local variables
    integer :: i,j,k

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. SLIP_WALL .or. bc(1,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-ng:lo(1)-1,j,k) = 0.5d0*(s(lo(1),j,k)+s(lo(1)-1,j,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. SLIP_WALL .or. bc(1,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1:hi(1)+ng,j,k) = 0.5d0*(s(hi(1),j,k)+s(hi(1)+1,j,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. SLIP_WALL .or. bc(2,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-ng:lo(2)-1,k) = 0.5d0*(s(i,lo(2),k)+s(i,lo(2)-1,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. SLIP_WALL .or. bc(2,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1:hi(2)+ng,k) = 0.5d0*(s(i,hi(2),k)+s(i,hi(2)+1,k))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. SLIP_WALL .or. bc(3,1) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-ng:lo(3)-1) = 0.5d0*(s(i,j,lo(3))+s(i,j,lo(3)-1))
          end do
       end do
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. SLIP_WALL .or. bc(3,2) .eq. NO_SLIP_WALL) then
       ! average interior and ghost into ghost
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1:hi(3)+ng) = 0.5d0*(s(i,j,hi(3))+s(i,j,hi(3)+1))
          end do
       end do
    end if

  end subroutine coefbc_3d

end module multifab_coefbc_module
