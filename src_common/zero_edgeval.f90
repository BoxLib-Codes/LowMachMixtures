module zero_edgeval_module

  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module

  implicit none

  private

  public :: zero_edgeval_physical, zero_edgeval_walls

contains

  subroutine zero_edgeval_physical(edge,start_comp,num_comp,the_bc_level)

    ! vel_bc_n(nlevs,dm) are the normal velocities

    type(multifab) , intent(inout) :: edge(:)
    integer        , intent(in   ) :: start_comp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(edge(1))),hi(get_dim(edge(1)))
    integer                  :: ng_e,i,dm,comp
    real(kind=dp_t), pointer :: epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:)

    ng_e = nghost(edge(1))
    dm = get_dim(edge(1))
    
    do i=1,nfabs(edge(1))
       epx => dataptr(edge(1),i)
       epy => dataptr(edge(2),i)
       lo = lwb(get_box(edge(1),i))
       hi = upb(get_box(edge(2),i))
       do comp=start_comp,start_comp+num_comp-1
          select case (dm)
          case (2)
             call zero_edgeval_physical_2d(epx(:,:,1,comp), epy(:,:,1,comp), ng_e, lo, hi, &
                                           the_bc_level%phys_bc_level_array(i,:,:))
          case (3)
             epz => dataptr(edge(3),i)
             call zero_edgeval_physical_3d(epx(:,:,:,comp), epy(:,:,:,comp), &
                                           epz(:,:,:,comp), ng_e, lo, hi, &
                                           the_bc_level%phys_bc_level_array(i,:,:))
          end select
       end do
    end do
 
  end subroutine zero_edgeval_physical

  subroutine zero_edgeval_physical_2d(edgex,edgey,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .ne. PERIODIC .and. bc(1,1) .ne. INTERIOR) then
       edgex(lo(1),lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .ne. PERIODIC .and. bc(1,2) .ne. INTERIOR) then
       edgex(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .ne. PERIODIC .and. bc(2,1) .ne. INTERIOR) then
       edgey(lo(1):hi(1),lo(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .ne. PERIODIC .and. bc(2,2) .ne. INTERIOR) then
       edgey(lo(1):hi(1),hi(2)+1) = 0.d0
    end if

  end subroutine zero_edgeval_physical_2d

  subroutine zero_edgeval_physical_3d(edgex,edgey,edgez,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .ne. PERIODIC .and. bc(1,1) .ne. INTERIOR) then
       edgex(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .ne. PERIODIC .and. bc(1,2) .ne. INTERIOR) then
       edgex(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .ne. PERIODIC .and. bc(2,1) .ne. INTERIOR) then
       edgey(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .ne. PERIODIC .and. bc(2,2) .ne. INTERIOR) then
       edgey(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .ne. PERIODIC .and. bc(3,1) .ne. INTERIOR) then
       edgez(lo(1):hi(1),lo(2):hi(2),lo(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .ne. PERIODIC .and. bc(3,2) .ne. INTERIOR) then
       edgez(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

  end subroutine zero_edgeval_physical_3d

  subroutine zero_edgeval_walls(edge,start_comp,num_comp,the_bc_level)

    ! vel_bc_n(nlevs,dm) are the normal velocities

    type(multifab) , intent(inout) :: edge(:)
    integer        , intent(in   ) :: start_comp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level

    ! Local
    integer                  :: lo(get_dim(edge(1))),hi(get_dim(edge(1)))
    integer                  :: ng_e,i,dm,comp
    real(kind=dp_t), pointer :: epx(:,:,:,:), epy(:,:,:,:), epz(:,:,:,:)

    ng_e = nghost(edge(1))
    dm = get_dim(edge(1))
    
    do i=1,nfabs(edge(1))
       epx => dataptr(edge(1),i)
       epy => dataptr(edge(2),i)
       lo = lwb(get_box(edge(1),i))
       hi = upb(get_box(edge(2),i))
       do comp=start_comp,start_comp+num_comp-1
          select case (dm)
          case (2)
             call zero_edgeval_walls_2d(epx(:,:,1,comp), epy(:,:,1,comp), ng_e, lo, hi, &
                                           the_bc_level%phys_bc_level_array(i,:,:))
          case (3)
             epz => dataptr(edge(3),i)
             call zero_edgeval_walls_3d(epx(:,:,:,comp), epy(:,:,:,comp), &
                                           epz(:,:,:,comp), ng_e, lo, hi, &
                                           the_bc_level%phys_bc_level_array(i,:,:))
          end select
       end do
    end do
 
  end subroutine zero_edgeval_walls

  subroutine zero_edgeval_walls_2d(edgex,edgey,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_WALL .or. bc(1,1) .eq. SLIP_WALL) then
       edgex(lo(1),lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_WALL .or. bc(1,2) .eq. SLIP_WALL) then
       edgex(hi(1)+1,lo(2):hi(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_WALL .or. bc(2,1) .eq. SLIP_WALL) then
       edgey(lo(1):hi(1),lo(2)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_WALL .or. bc(2,2) .eq. SLIP_WALL) then
       edgey(lo(1):hi(1),hi(2)+1) = 0.d0
    end if

  end subroutine zero_edgeval_walls_2d

  subroutine zero_edgeval_walls_3d(edgex,edgey,edgez,ng_e,lo,hi,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_e
    real(kind=dp_t), intent(inout) :: edgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(inout) :: edgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    integer        , intent(in   ) :: bc(:,:)

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. NO_SLIP_WALL .or. bc(1,1) .eq. SLIP_WALL) then
       edgex(lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. NO_SLIP_WALL .or. bc(1,2) .eq. SLIP_WALL) then
       edgex(hi(1)+1,lo(2):hi(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. NO_SLIP_WALL .or. bc(2,1) .eq. SLIP_WALL) then
       edgey(lo(1):hi(1),lo(2),lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. NO_SLIP_WALL .or. bc(2,2) .eq. SLIP_WALL) then
       edgey(lo(1):hi(1),hi(2)+1,lo(3):hi(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. NO_SLIP_WALL .or. bc(3,1) .eq. SLIP_WALL) then
       edgez(lo(1):hi(1),lo(2):hi(2),lo(3)) = 0.d0
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. NO_SLIP_WALL .or. bc(3,2) .eq. SLIP_WALL) then
       edgez(lo(1):hi(1),lo(2):hi(2),hi(3)+1) = 0.d0
    end if

  end subroutine zero_edgeval_walls_3d

end module zero_edgeval_module
