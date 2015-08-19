module multifab_physbc_extrap_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module

  implicit none

  private

  public :: multifab_physbc_extrap

contains


  ! multifab_physbc:        Fills in all ghost cells to the same value, which is the value AT the boundary.
  !                         HOEXTRAP/FOEXTRAP uses boundary conditions (Neumann) and 1/2 interior points.
  !                         EXT_DIR copies the supplied Dirichlet condition into the ghost cells.

  ! multifab_physbc_extrap: Fills in two cell-centered ghost cells by using cubic extrapolation using the
  !                         boundary condition and three interior points.  We are computing the values 
  !                         extrapolated to the ghost cell-centers.  For FOEXTRAP/HOEXTRAP for now
  !                         we assume homogeneous Neumann.  For EXT_DIR, we assume the ghost cell on input
  !                         contains the Dirichlet value.

  subroutine multifab_physbc_extrap(s,start_scomp,start_bccomp,num_comp,the_bc_level, &
                                    time_in,dx_in,prob_lo_in,prob_hi_in,increment_bccomp_in)

    ! this fills 2 ghost cells for conetration using cubic interpolation
    ! here, we extrapolate to the ghost cell-centers, not the physical boundary location
    ! for Dirichlet conditions, the ghost cell on input contains the Dirichlet value
    ! for Neumann conditions, we assume homogeneous

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: start_scomp,start_bccomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ), optional :: time_in,dx_in(:),prob_lo_in(:),prob_hi_in(:)
    logical        ,  intent(in  ), optional :: increment_bccomp_in
   
    ! Local
    integer :: lo(get_dim(s)),hi(get_dim(s))
    integer :: i,ng,dm,scomp,bccomp
    logical :: increment_bccomp

    real(kind=dp_t), allocatable :: dx(:)
    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"multifab_physbc_extrap")

    ng = nghost(s)
    dm = get_dim(s)

    allocate(dx(dm))

    if (present(dx_in)) then
       dx(1:dm) = dx_in(1:dm)
    else
       dx = 1.d0
    end if

    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       do scomp=start_scomp,start_scomp+num_comp-1
          if (increment_bccomp) then
             bccomp = start_bccomp + (scomp-start_scomp)
          else
             bccomp = start_bccomp
          end if
          select case (dm)
          case (2)
             call physbc_extrap_2d(sp(:,:,1,scomp), lo, hi, ng, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          case (3)
             call physbc_extrap_3d(sp(:,:,:,scomp), lo, hi, ng, &
                                the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          end select
       end do
    end do

    deallocate(dx)
 
    call destroy(bpt)

  end subroutine multifab_physbc_extrap

  subroutine physbc_extrap_2d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j
    real(kind=dp_t) :: fac,dir_val

    fac = 1.d0/29.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       ! (NOTE: The full stencils with the Neumann value are:
       !  fac*(- 39*Neumann_val*dx + 15...
       !  fac*(-264*Neumann_val*dx + 15...
       ! The sign in front of the Neumann_val changes on the lo/hi side
       do j=lo(2)-ng,hi(2)+ng
          s(lo(1)-1,j) = fac*(  15.d0*s(lo(1),j) +  18.d0*s(lo(1)+1,j) -  4.d0*s(lo(1)+2,j))
          s(lo(1)-2,j) = fac*(-242.d0*s(lo(1),j) + 336.d0*s(lo(1)+1,j) - 65.d0*s(lo(1)+2,j))
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do j=lo(2)-ng,hi(2)+ng
          dir_val = s(lo(1)-1,j)
          s(lo(1)-1,j) =  -26.d0*dir_val +  44.d0*s(lo(1),j) -  20.d0*s(lo(1)+1,j) +  3.d0*s(lo(1)+2,j)
          s(lo(1)-2,j) = -176.d0*dir_val + 286.d0*s(lo(1),j) - 128.d0*s(lo(1)+1,j) + 19.d0*s(lo(1)+2,j)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do j=lo(2)-ng,hi(2)+ng
          s(hi(1)+1,j) = fac*(  15.d0*s(hi(1),j) +  18.d0*s(hi(1)-1,j) -  4.d0*s(hi(1)-2,j))
          s(hi(1)+2,j) = fac*(-242.d0*s(hi(1),j) + 336.d0*s(hi(1)-1,j) - 65.d0*s(hi(1)-2,j))
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do j=lo(2)-ng,hi(2)+ng
          dir_val = s(hi(1)+1,j)
          s(hi(1)+1,j) =  -26.d0*dir_val +  44.d0*s(hi(1),j) -  20.d0*s(hi(1)-1,j) +  3.d0*s(hi(1)-2,j)
          s(hi(1)+2,j) = -176.d0*dir_val + 286.d0*s(hi(1),j) - 128.d0*s(hi(1)-1,j) + 19.d0*s(hi(1)-2,j)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do i=lo(1)-ng,hi(1)+ng
          s(i,lo(2)-1) = fac*(  15.d0*s(i,lo(2)) +  18.d0*s(i,lo(2)+1) -  4.d0*s(i,lo(2)+2))
          s(i,lo(2)-2) = fac*(-242.d0*s(i,lo(2)) + 336.d0*s(i,lo(2)+1) - 65.d0*s(i,lo(2)+2))
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do i=lo(1)-ng,hi(1)+ng
          dir_val = s(i,lo(2)-1)
          s(i,lo(2)-1) =  -26.d0*dir_val +  44.d0*s(i,lo(2)) -  20.d0*s(i,lo(2)+1) +  3.d0*s(i,lo(2)+2)
          s(i,lo(2)-2) = -176.d0*dir_val + 286.d0*s(i,lo(2)) - 128.d0*s(i,lo(2)+1) + 19.d0*s(i,lo(2)+2)
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do i=lo(1)-ng,hi(1)+ng
          s(i,hi(2)+1) = fac*(  15.d0*s(i,hi(2)) +  18.d0*s(i,hi(2)-1) -  4.d0*s(i,hi(2)-2))
          s(i,hi(2)+2) = fac*(-242.d0*s(i,hi(2)) + 336.d0*s(i,hi(2)-1) - 65.d0*s(i,hi(2)-2))
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do i=lo(1)-ng,hi(1)+ng
          dir_val = s(i,hi(2)+1)
          s(i,hi(2)+1) =  -26.d0*dir_val +  44.d0*s(i,hi(2)) -  20.d0*s(i,hi(2)-1) +  3.d0*s(i,hi(2)-2)
          s(i,hi(2)+2) = -176.d0*dir_val + 286.d0*s(i,hi(2)) - 128.d0*s(i,hi(2)-1) + 19.d0*s(i,hi(2)-2)
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_extrap_2d

  subroutine physbc_extrap_3d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k
    real(kind=dp_t) :: fac,dir_val

    fac = 1.d0/29.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(lo(1)-1,j,k) = fac*(  15.d0*s(lo(1),j,k) +  18.d0*s(lo(1)+1,j,k) -  4.d0*s(lo(1)+2,j,k))
             s(lo(1)-2,j,k) = fac*(-242.d0*s(lo(1),j,k) + 336.d0*s(lo(1)+1,j,k) - 65.d0*s(lo(1)+2,j,k))
          end do
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             dir_val = s(lo(1)-1,j,k)
             s(lo(1)-1,j,k) =  -26.d0*dir_val +  44.d0*s(lo(1),j,k) -  20.d0*s(lo(1)+1,j,k) +  3.d0*s(lo(1)+2,j,k)
             s(lo(1)-2,j,k) = -176.d0*dir_val + 286.d0*s(lo(1),j,k) - 128.d0*s(lo(1)+1,j,k) + 19.d0*s(lo(1)+2,j,k)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             s(hi(1)+1,j,k) = fac*(  15.d0*s(hi(1),j,k) +  18.d0*s(hi(1)-1,j,k) -  4.d0*s(hi(1)-2,j,k))
             s(hi(1)+2,j,k) = fac*(-242.d0*s(hi(1),j,k) + 336.d0*s(hi(1)-1,j,k) - 65.d0*s(hi(1)-2,j,k))
          end do
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do k=lo(3)-ng,hi(3)+ng
          do j=lo(2)-ng,hi(2)+ng
             dir_val = s(hi(1)+1,j,k)
             s(hi(1)+1,j,k) =  -26.d0*dir_val +  44.d0*s(hi(1),j,k) -  20.d0*s(hi(1)-1,j,k) +  3.d0*s(hi(1)-2,j,k)
             s(hi(1)+2,j,k) = -176.d0*dir_val + 286.d0*s(hi(1),j,k) - 128.d0*s(hi(1)-1,j,k) + 19.d0*s(hi(1)-2,j,k)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,lo(2)-1,k) = fac*(  15.d0*s(i,lo(2),k) +  18.d0*s(i,lo(2)+1,k) -  4.d0*s(i,lo(2)+2,k))
             s(i,lo(2)-2,k) = fac*(-242.d0*s(i,lo(2),k) + 336.d0*s(i,lo(2)+1,k) - 65.d0*s(i,lo(2)+2,k))
          end do
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             dir_val = s(i,lo(2)-1,k)
             s(i,lo(2)-1,k) =  -26.d0*dir_val +  44.d0*s(i,lo(2),k) -  20.d0*s(i,lo(2)+1,k) +  3.d0*s(i,lo(2)+2,k)
             s(i,lo(2)-2,k) = -176.d0*dir_val + 286.d0*s(i,lo(2),k) - 128.d0*s(i,lo(2)+1,k) + 19.d0*s(i,lo(2)+2,k)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,hi(2)+1,k) = fac*(  15.d0*s(i,hi(2),k) +  18.d0*s(i,hi(2)-1,k) -  4.d0*s(i,hi(2)-2,k))
             s(i,hi(2)+2,k) = fac*(-242.d0*s(i,hi(2),k) + 336.d0*s(i,hi(2)-1,k) - 65.d0*s(i,hi(2)-2,k))
          end do
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do k=lo(3)-ng,hi(3)+ng
          do i=lo(1)-ng,hi(1)+ng
             dir_val = s(i,hi(2)+1,k)
             s(i,hi(2)+1,k) =  -26.d0*dir_val +  44.d0*s(i,hi(2),k) -  20.d0*s(i,hi(2)-1,k) +  3.d0*s(i,hi(2)-2,k)
             s(i,hi(2)+2,k) = -176.d0*dir_val + 286.d0*s(i,hi(2),k) - 128.d0*s(i,hi(2)-1,k) + 19.d0*s(i,hi(2)-2,k)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,lo(3)-1) = fac*(  15.d0*s(i,j,lo(3)) +  18.d0*s(i,j,lo(3)+1) -  4.d0*s(i,j,lo(3)+2))
             s(i,j,lo(3)-2) = fac*(-242.d0*s(i,j,lo(3)) + 336.d0*s(i,j,lo(3)+1) - 65.d0*s(i,j,lo(3)+2))
          end do
       end do
    else if (bc(3,1) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             dir_val = s(i,j,lo(3)-1)
             s(i,j,lo(3)-1) =  -26.d0*dir_val +  44.d0*s(i,j,lo(3)) -  20.d0*s(i,j,lo(3)+1) +  3.d0*s(i,j,lo(3)+2)
             s(i,j,lo(3)-2) = -176.d0*dir_val + 286.d0*s(i,j,lo(3)) - 128.d0*s(i,j,lo(3)+1) + 19.d0*s(i,j,lo(3)+2)
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP) then
       ! cubic extrapolation using 3 interior cell averages and Neumann bc
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             s(i,j,hi(3)+1) = fac*(  15.d0*s(i,j,hi(3)) +  18.d0*s(i,j,hi(3)-1) -  4.d0*s(i,j,hi(3)-2))
             s(i,j,hi(3)+2) = fac*(-242.d0*s(i,j,hi(3)) + 336.d0*s(i,j,hi(3)-1) - 65.d0*s(i,j,hi(3)-2))
          end do
       end do
    else if (bc(3,2) .eq. EXT_DIR) then
       ! cubic extrapolation using 3 interior cell averages and Dirichlet bc
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             dir_val = s(i,j,hi(3)+1)
             s(i,j,hi(3)+1) =  -26.d0*dir_val +  44.d0*s(i,j,hi(3)) -  20.d0*s(i,j,hi(3)-1) +  3.d0*s(i,j,hi(3)-2)
             s(i,j,hi(3)+2) = -176.d0*dir_val + 286.d0*s(i,j,hi(3)) - 128.d0*s(i,j,hi(3)-1) + 19.d0*s(i,j,hi(3)-2)
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_extrap_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_extrap_3d

end module multifab_physbc_extrap_module
