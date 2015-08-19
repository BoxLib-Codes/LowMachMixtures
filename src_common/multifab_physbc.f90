module multifab_physbc_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: multifab_physbc

contains

  ! multifab_physbc:        Fills in all ghost cells to the same value, which is the value AT the boundary.
  !                         HOEXTRAP/FOEXTRAP uses boundary conditions (Neumann) and 1/2 interior points.
  !                         EXT_DIR copies the supplied Dirichlet condition into the ghost cells.

  ! multifab_physbc_extrap: Fills in two cell-centered ghost cells by using cubic extrapolation using the
  !                         boundary condition and three interior points.  We are computing the values 
  !                         extrapolated to the ghost cell-centers.  For FOEXTRAP/HOEXTRAP for now
  !                         we assume homogeneous Neumann.  For EXT_DIR, we assume the ghost cell on input
  !                         contains the Dirichlet value.

  subroutine multifab_physbc(s,start_scomp,start_bccomp,num_comp,the_bc_level, &
                             time_in,dx_in,prob_lo_in,prob_hi_in,increment_bccomp_in)

    ! this fills ghost cells for rho and pressure/phi.
    ! as well as for transport coefficients (alpha/beta/gamma)

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

    call build(bpt,"multifab_physbc")

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
             call physbc_2d(sp(:,:,1,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          case (3)
             call physbc_3d(sp(:,:,:,scomp), lo, hi, ng, &
                            the_bc_level%adv_bc_level_array(i,:,:,bccomp),bccomp, dx)
          end select
       end do
    end do

    deallocate(dx)
 
    call destroy(bpt)

  end subroutine multifab_physbc

  subroutine physbc_2d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j
    real(kind=dp_t) :: x,y,threeEighths,nineEighths,oneEighth

    threeEighths = 3.d0/8.d0
    nineEighths = 9.d0/8.d0
    oneEighth = 1.d0/8.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       x = prob_lo(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(lo(1)-ng:lo(1)-1,j) = s(lo(1),j) - 0.5d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       x = prob_lo(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(lo(1)-ng:lo(1)-1,j) = -(threeEighths)*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
               + (nineEighths)*s(lo(1),j) - (oneEighth)*s(lo(1)+1,j)
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_lo(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(lo(1)-ng:lo(1)-1,j) = inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       x = prob_hi(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(hi(1)+1:hi(1)+ng,j) = s(hi(1),j) + 0.5d0*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       x = prob_hi(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(hi(1)+1:hi(1)+ng,j) = (threeEighths)*dx(1)*inhomogeneous_bc_val_2d(bccomp,x,y) &
               + (nineEighths)*s(hi(1),j) - (oneEighth)*s(hi(1)-1,j)
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_hi(1)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          s(hi(1)+1:hi(1)+ng,j) = inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       y = prob_lo(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,lo(2)-ng:lo(2)-1) = s(i,lo(2)) - 0.5d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       y = prob_lo(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,lo(2)-ng:lo(2)-1) = -(threeEighths)*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
               + (nineEighths)*s(i,lo(2)) - (oneEighth)*s(i,lo(2)+1)
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       y = prob_lo(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,lo(2)-ng:lo(2)-1) = inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       y = prob_hi(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,hi(2)+1:hi(2)+ng) = s(i,hi(2)) + 0.5d0*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       y = prob_hi(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,hi(2)+1:hi(2)+ng) = (threeEighths)*dx(2)*inhomogeneous_bc_val_2d(bccomp,x,y) &
               + (nineEighths)*s(i,hi(2)) - (oneEighth)*s(i,hi(2)-1)
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       y = prob_hi(2)
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          s(i,hi(2)+1:hi(2)+ng) = inhomogeneous_bc_val_2d(bccomp,x,y)
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_2d

  subroutine physbc_3d(s,lo,hi,ng,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z,threeEighths,nineEighths,oneEighth

    threeEighths = 3.d0/8.d0
    nineEighths = 9.d0/8.d0
    oneEighth = 1.d0/8.d0

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       x = prob_lo(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1)-ng:lo(1)-1,j,k) = s(lo(1),j,k) &
                  - 0.5d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(1,1) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       x = prob_lo(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1)-ng:lo(1)-1,j,k) = &
                  -(threeEighths)*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(lo(1),j,k) - (oneEighth)*s(lo(1)+1,j,k)
          end do
       end do
    else if (bc(1,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_lo(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(lo(1)-ng:lo(1)-1,j,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       x = prob_hi(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1:hi(1)+ng,j,k) = s(hi(1),j,k) &
                  + 0.5d0*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(1,2) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       x = prob_hi(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1:hi(1)+ng,j,k) = &
                  (threeEighths)*dx(1)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(hi(1),j,k) - (oneEighth)*s(hi(1)-1,j,k)
          end do
       end do
    else if (bc(1,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       x = prob_hi(1)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2)-ng,hi(2)+ng
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             s(hi(1)+1:hi(1)+ng,j,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       y = prob_lo(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)-ng:lo(2)-1,k) = s(i,lo(2),k) &
                  - 0.5d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(2,1) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       y = prob_lo(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)-ng:lo(2)-1,k) = &
                  -(threeEighths)*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(i,lo(2),k) - (oneEighth)*s(i,lo(2)+1,k)
          end do
       end do
    else if (bc(2,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       y = prob_lo(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,lo(2)-ng:lo(2)-1,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       y = prob_hi(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1:hi(2)+ng,k) = s(i,hi(2),k) &
                  + 0.5d0*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(2,2) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       y = prob_hi(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1:hi(2)+ng,k) = &
                  (threeEighths)*dx(2)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(i,hi(2),k) - (oneEighth)*s(i,hi(2)-1,k)
          end do
       end do
    else if (bc(2,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       y = prob_hi(2)
       do k=lo(3)-ng,hi(3)+ng
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,hi(2)+1:hi(2)+ng,k) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       z = prob_lo(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,lo(3)-ng:lo(3)-1) = s(i,j,lo(3)) &
                  - 0.5d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(3,1) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       z = prob_lo(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,lo(3)-ng:lo(3)-1) = &
                  -(threeEighths)*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(i,j,lo(3)) - (oneEighth)*s(i,j,lo(3)+1)
          end do
       end do
    else if (bc(3,1) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       z = prob_lo(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,lo(3)-ng:lo(3)-1) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. FOEXTRAP) then
       ! linear extrapolation using interior point and Neumann bc
       z = prob_hi(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,hi(3)+1:hi(3)+ng) = s(i,j,hi(3)) &
                  + 0.5d0*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(3,2) .eq. HOEXTRAP) then
       ! quadratric extrapolation using two interior points and Neumann bc
       z = prob_hi(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,hi(3)+1:hi(3)+ng) = &
                  (threeEighths)*dx(3)*inhomogeneous_bc_val_3d(bccomp,x,y,z) &
                  + (nineEighths)*s(i,j,hi(3)) - (oneEighth)*s(i,j,hi(3)-1)
          end do
       end do
    else if (bc(3,2) .eq. EXT_DIR) then
       ! dirichlet condition obtained from inhomogeneous_bc_val
       z = prob_hi(3)
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             s(i,j,hi(3)+1:hi(3)+ng) = inhomogeneous_bc_val_3d(bccomp,x,y,z)
          end do
       end do
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_3d

end module multifab_physbc_module
