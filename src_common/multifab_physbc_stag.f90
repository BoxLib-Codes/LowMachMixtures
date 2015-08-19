module multifab_physbc_stag_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use bc_module
  use bl_error_module
  use inhomogeneous_bc_val_module
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

  private

  public :: multifab_physbc_domainvel, multifab_physbc_macvel, &
       set_inhomogeneous_vel_bcs, modify_traction_bcs

contains

  subroutine multifab_physbc_domainvel(s,bccomp,the_bc_level,dx,vel_bc_n)

    ! vel_bc_n(nlevs,dm) are the normal velocities

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: bccomp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(multifab) , intent(in   ), optional :: vel_bc_n(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,ng_v,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:), vp(:,:,:,:)
    logical                  :: use_inhomogeneous

    type(bl_prof_timer),save :: bpt

    call build(bpt,"multifab_physbc_domainvel")

    use_inhomogeneous = .false.
    if (present(vel_bc_n)) then
       use_inhomogeneous = .true.
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_domainvel expects bccomp <= dm')
    end if

    ng_s = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          if (use_inhomogeneous) then
             ng_v = vel_bc_n(bccomp)%ng
             vp => dataptr(vel_bc_n(bccomp),i)
             call physbc_domainvel_2d_inhomogeneous(sp(:,:,1,1), ng_s, &
                                                    vp(:,:,1,1), ng_v, lo, hi, &
                                                    the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                    bccomp,dx)
          else
             call physbc_domainvel_2d(sp(:,:,1,1), ng_s, lo, hi, &
                                      the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                      bccomp,dx)
          end if
       case (3)
          if (use_inhomogeneous) then
             ng_v = vel_bc_n(bccomp)%ng
             vp => dataptr(vel_bc_n(bccomp),i)
             call physbc_domainvel_3d_inhomogeneous(sp(:,:,:,1), ng_s, &
                                                    vp(:,:,:,1), ng_v, lo, hi, &
                                                    the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                    bccomp,dx)
          else
             call physbc_domainvel_3d(sp(:,:,:,1), ng_s, lo, hi, &
                                      the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                      bccomp,dx)
          end if
       end select
    end do
 
    call destroy(bpt)

  end subroutine multifab_physbc_domainvel

  subroutine physbc_domainvel_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_domainvel_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          ! need to set ghost cells behind to avoid propagating NaNs in the bds scheme
          ! that should be wiped out (but don't since zero times a NaN is a Nan)
          s(lo(1)-ng_s:lo(1),lo(2):hi(2)) = 0.d0
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(hi(1)+1:hi(1)+1+ng_s,lo(2):hi(2)) = 0.d0
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1):hi(1),lo(2)-ng_s:lo(2)) = 0.d0
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1):hi(1),hi(2)+1:hi(2)+1+ng_s) = 0.d0
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

  end subroutine physbc_domainvel_2d

  subroutine physbc_domainvel_2d_inhomogeneous(s,ng_s,vel_bc_n,ng_v,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_v
    real(kind=dp_t), intent(inout) ::      s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(inout) :: vel_bc_n(lo(1)-ng_v:,lo(2)-ng_v:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_domainvel_2d_inhomogeneous requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          ! need to set ghost cells behind to avoid propagating NaNs in the bds scheme
          ! that should be wiped out (but don't since zero times a NaN is a Nan)
          do j=lo(2),hi(2)
             s(lo(1)-ng_s:lo(1),j) = vel_bc_n(lo(1),j)
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d_inhomogeneous: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do j=lo(2),hi(2)
             s(hi(1)+1:hi(1)+1+ng_s,j) = vel_bc_n(hi(1)+1,j)
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d_inhomogeneous: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do i=lo(1),hi(1)
             s(i,lo(2)-ng_s:lo(2)) = vel_bc_n(i,lo(2))
          end do
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d_inhomogeneous: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do i=lo(1),hi(1)
             s(i,hi(2)+1:hi(2)+1+ng_s) = vel_bc_n(i,hi(2)+1)
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_2d_inhomogeneous: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

  end subroutine physbc_domainvel_2d_inhomogeneous

  subroutine physbc_domainvel_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_domainvel_3d requires bccomp = 1, 2, or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          ! need to set ghost cells behind to avoid propagating NaNs in the bds scheme
          ! that should be wiped out (but don't since zero times a NaN is a Nan)
          s(lo(1)-ng_s:lo(1),lo(2):hi(2),lo(3):hi(3)) = 0.d0
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(hi(1)+1:hi(1)+1+ng_s,lo(2):hi(2),lo(3):hi(3)) = 0.d0
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1):hi(1),lo(2)-ng_s:lo(2),lo(3):hi(3)) = 0.d0
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1):hi(1),hi(2)+1:hi(2)+1+ng_s,lo(3):hi(3)) = 0.d0
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    ! staggered z-velocity
    if (bccomp .eq. 3) then
       if (bc(3,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          s(lo(1):hi(1),lo(2):hi(2),lo(3)-ng_s:lo(3)) = 0.d0
       else if (bc(3,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
      end if
   end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

   ! staggered z-velocity
   if (bccomp .eq. 3) then
      if (bc(3,2) .eq. DIR_VEL) then
         ! set domain face value to Dirichlet value
         s(lo(1):hi(1),lo(2):hi(2),hi(3)+1:hi(3)+1+ng_s) = 0.d0
      else if (bc(3,2) .eq. INTERIOR) then
         ! either periodic or interior; do nothing
      else
         print *,'physbc_domainvel_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
         call bl_error('NOT SUPPORTED')
      end if
   end if

 end subroutine physbc_domainvel_3d

  subroutine physbc_domainvel_3d_inhomogeneous(s,ng_s,vel_bc_n,ng_v,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_v
    real(kind=dp_t), intent(inout) ::      s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: vel_bc_n(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! local
    integer :: i,j,k

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_domainvel_3d_inhomogeneous requires bccomp = 1, 2, or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          ! need to set ghost cells behind to avoid propagating NaNs in the bds scheme
          ! that should be wiped out (but don't since zero times a NaN is a Nan)
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                s(lo(1)-ng_s:lo(1),j,k) = vel_bc_n(lo(1),j,k)
             end do
          end do
       else if (bc(1,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d_inhomogeneous: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    ! staggered x-velocity
    if (bccomp .eq. 1) then
       if (bc(1,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                s(hi(1)+1:hi(1)+1+ng_s,j,k) = vel_bc_n(hi(1)+1,j,k)
             end do
          end do
       else if (bc(1,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d_inhomogeneous: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)
                s(i,lo(2)-ng_s:lo(2),k) = vel_bc_n(i,lo(2),k)
             end do
          end do
       else if (bc(2,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d_inhomogeneous: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    ! staggered y-velocity
    if (bccomp .eq. 2) then
       if (bc(2,2) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)
                s(i,hi(2)+1:hi(2)+1+ng_s,k) = vel_bc_n(i,hi(2)+1,k)
             end do
          end do
       else if (bc(2,2) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d_inhomogeneous: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    ! staggered z-velocity
    if (bccomp .eq. 3) then
       if (bc(3,1) .eq. DIR_VEL) then
          ! set domain face value to Dirichlet value
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s(i,j,lo(3)-ng_s:lo(3)) = vel_bc_n(i,j,lo(3))
             end do
          end do
       else if (bc(3,1) .eq. INTERIOR) then
          ! either periodic or interior; do nothing
       else
          print *,'physbc_domainvel_3d_inhomogeneous: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
      end if
   end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

   ! staggered z-velocity
   if (bccomp .eq. 3) then
      if (bc(3,2) .eq. DIR_VEL) then
         ! set domain face value to Dirichlet value
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                s(i,j,hi(3)+1:hi(3)+1+ng_s) = vel_bc_n(i,j,hi(3)+1)
             end do
          end do
      else if (bc(3,2) .eq. INTERIOR) then
         ! either periodic or interior; do nothing
      else
         print *,'physbc_domainvel_3d_inhomogeneous: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
         call bl_error('NOT SUPPORTED')
      end if
   end if

 end subroutine physbc_domainvel_3d_inhomogeneous

  subroutine multifab_physbc_macvel(s,bccomp,the_bc_level,dx,vel_bc_t)

    ! in 2D, vel_bc_t(nlevs,2) respresents
    !   1. y-velocity bc on x-faces (nodal)
    !   2. x-velocity bc on y-faces (nodal)
    ! in 3D, vel_bc_t(nlevs,6) represents
    !   1. y-velocity bc on x-faces (nodal in y and x)
    !   2. z-velocity bc on x-faces (nodal in z and x)
    !   3. x-velocity bc on y-faces (nodal in x and y)
    !   4. z-velocity bc on y-faces (nodal in z and y)
    !   5. x-velocity bc on z-faces (nodal in x and z)
    !   6. y-velocity bc on z-faces (nodal in y and z)

    type(multifab) , intent(inout) :: s
    integer        , intent(in   ) :: bccomp
    type(bc_level) , intent(in   ) :: the_bc_level
    real(kind=dp_t), intent(in   ) :: dx(:)
    type(multifab) , intent(in   ), optional :: vel_bc_t(:)

    ! Local
    integer                  :: lo(get_dim(s)),hi(get_dim(s))
    integer                  :: i,ng_s,ng_v,dm
    real(kind=dp_t), pointer :: sp(:,:,:,:), vp1(:,:,:,:), vp2(:,:,:,:)
    logical                  :: use_inhomogeneous

    type(bl_prof_timer),save :: bpt

    call build(bpt,"multifab_physbc_macvel")

    use_inhomogeneous = .false.
    if (present(vel_bc_t)) then
       use_inhomogeneous = .true.
    end if

    if (bccomp .gt. get_dim(s)) then
       call bl_error('multifab_physbc_macvel expects bccomp <= dm')
    end if

    ng_s = nghost(s)
    dm = get_dim(s)
    
    do i=1,nfabs(s)
       sp => dataptr(s,i)
       lo = lwb(get_box(s,i))
       hi = upb(get_box(s,i))
       select case (dm)
       case (2)
          if (use_inhomogeneous) then
             ng_v = vel_bc_t(bccomp)%ng
             if (bccomp .eq. 1) then
                vp1 => dataptr(vel_bc_t(2),i)
             else if (bccomp .eq. 2) then
                vp1 => dataptr(vel_bc_t(1),i)
             end if
             call physbc_macvel_2d_inhomogeneous(sp(:,:,1,1), ng_s, &
                                                 vp1(:,:,1,1), ng_v, lo, hi, &
                                                 the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                 bccomp,dx)
          else
             call physbc_macvel_2d(sp(:,:,1,1), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                   bccomp,dx)
          end if
       case (3)
          if (use_inhomogeneous) then
             ng_v = vel_bc_t(bccomp)%ng
             if (bccomp .eq. 1) then
                vp1 => dataptr(vel_bc_t(3),i)
                vp2 => dataptr(vel_bc_t(5),i)
             else if (bccomp .eq. 2) then
                vp1 => dataptr(vel_bc_t(1),i)
                vp2 => dataptr(vel_bc_t(6),i)
             else if (bccomp .eq. 3) then
                vp1 => dataptr(vel_bc_t(2),i)
                vp2 => dataptr(vel_bc_t(4),i)
             end if
             call physbc_macvel_3d_inhomogeneous(sp(:,:,:,1), ng_s, &
                                                 vp1(:,:,:,1), vp2(:,:,:,1), ng_v, lo, hi, &
                                                 the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                                 bccomp,dx)
          else
             call physbc_macvel_3d(sp(:,:,:,1), ng_s, lo, hi, &
                                   the_bc_level%adv_bc_level_array(i,:,:,bccomp), &
                                   bccomp,dx)
          end if
       end select
    end do
 
    call destroy(bpt)

  end subroutine multifab_physbc_macvel

  subroutine physbc_macvel_2d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_macvel_2d requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = -s(lo(1),j)
             ! higher-order stencil
             ! s(lo(1)-ng_s:lo(1)-1,j) = -2.d0*s(lo(1),j) + (1.d0/3.d0)*s(lo(1)+1,j)
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = s(lo(1),j)
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = -s(hi(1),j)
             ! higher-order stencil
             ! s(hi(1)+1:hi(1)+ng_s,j) = -2.d0*s(hi(1),j) + (1.d0/3.d0)*s(hi(1)-1,j)
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = s(hi(1),j)
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = -s(i,lo(2))
             ! higher-order stencil
             ! s(i,lo(2)-ng_s:lo(2)-1) = -2.d0*s(i,lo(2)) + (1.d0/3.d0)*s(i,lo(2)+1)
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = s(i,lo(2))
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = -s(i,hi(2))
             ! higher-order stencil
             ! s(i,hi(2)+1:hi(2)+ng_s) = -2.d0*s(i,hi(2)) + (1.d0/3.d0)*s(i,hi(2)-1)
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = s(i,hi(2))
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_2d

  subroutine physbc_macvel_2d_inhomogeneous(s,ng_s,v_bc,ng_v,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_v
    real(kind=dp_t), intent(inout) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(in   ) :: v_bc(lo(1)-ng_v:,lo(2)-ng_v:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j

    if (bccomp .ne. 1 .and. bccomp .ne. 2) then
       call bl_error('physbc_macvel_2d_inhomogeneous requires bccomp = 1 or 2')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = 2.d0*v_bc(lo(1),j) - s(lo(1),j)
             ! higher-order stencil
             ! s(lo(1)-ng_s:lo(1)-1,j) = (8.d0/3.d0)*v_bc(lo(1),j) - 2.d0*s(lo(1),j) + (1.d0/3.d0)*s(lo(1)+1,j)
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d_inhomogeneous: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)+1
             s(lo(1)-ng_s:lo(1)-1,j) = s(lo(1),j) - dx(1)*v_bc(lo(1),j)
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d_inhomogeneous: bc(1,1) =',bc(1,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = 2.d0*v_bc(hi(1)+1,j) - s(hi(1),j)
             ! higher-order stencil
             ! s(hi(1)+1:hi(1)+ng_s,j) = (8.d0/3.d0)*v_bc(hi(1)+1,j) - 2.d0*s(hi(1),j) + (1.d0/3.d0)*s(hi(1)-1,j)
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_2d_inhomogeneous: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)+1
             s(hi(1)+1:hi(1)+ng_s,j) = s(hi(1),j) + dx(1)*v_bc(hi(1)+1,j)
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d_inhomogeneous: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = 2.d0*v_bc(i,lo(2)) - s(i,lo(2))
             ! higher-order stencil
             ! s(i,lo(2)-ng_s:lo(2)-1) = (8.d0/3.d0)*v_bc(i,lo(2)) - 2.d0*s(i,lo(2)) + (1.d0/3.d0)*s(i,lo(2)+1)
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do i=lo(1),hi(1)+1
             s(i,lo(2)-ng_s:lo(2)-1) = s(i,lo(2)) - dx(1)*v_bc(i,lo(2))
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d_inhomogeneous: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d_inhomogeneous: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = 2.d0*v_bc(i,hi(2)+1) - s(i,hi(2))
             ! higher-order stencil
             ! s(i,hi(2)+1:hi(2)+ng_s) = (8.d0/3.d0)*v_bc(i,hi(2)+1) - 2.d0*s(i,hi(2)) + (1.d0/3.d0)*s(i,hi(2)-1)
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do i=lo(1),hi(1)+1
             s(i,hi(2)+1:hi(2)+ng_s) = s(i,hi(2)) + dx(1)*v_bc(i,hi(2)+1)
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_2d_inhomogeneous: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_2d_inhomogeneous: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_2d_inhomogeneous

  subroutine physbc_macvel_3d(s,ng_s,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s
    real(kind=dp_t), intent(inout) :: s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! Local variables
    integer :: i,j,k

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_macvel_3d requires bccomp = 1, 2 or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = -s(lo(1),j,k)
                ! higher-order stencil
                ! s(lo(1)-ng_s:lo(1)-1,j,k) = -2.d0*s(lo(1),j,k) + (1.d0/3.d0)*s(lo(1)+1,j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(lo(1)-ng_s:lo(1)-1,j,k) = -s(lo(1),j,k)
                ! higher-order stencil
                ! s(lo(1)-ng_s:lo(1)-1,j,k) = -2.d0*s(lo(1),j,k) + (1.d0/3.d0)*s(lo(1)+1,j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(1,1) =',bc(1,1),'for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = -s(hi(1),j,k)
                ! higher-order stencil
                ! s(hi(1)+1:hi(1)+ng_s,j,k) = -2.d0*s(hi(1),j,k) + (1.d0/3.d0)*s(hi(1)-1,j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(hi(1)+1:hi(1)+ng_s,j,k) = -s(hi(1),j,k)
                ! higher-order stencil
                ! s(hi(1)+1:hi(1)+ng_s,j,k) = -2.d0*s(hi(1),j,k) + (1.d0/3.d0)*s(hi(1)-1,j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = -s(i,lo(2),k)
                ! higher-order stencil
                ! s(i,lo(2)-ng_s:lo(2)-1,k) = -2.d0*s(i,lo(2),k) + (1.d0/3.d0)*s(i,lo(2)+1,k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,lo(2)-ng_s:lo(2)-1,k) = -s(i,lo(2),k)
                ! higher-order stencil
                ! s(i,lo(2)-ng_s:lo(2)-1,k) = -2.d0*s(i,lo(2),k) + (1.d0/3.d0)*s(i,lo(2)+1,k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = -s(i,hi(2),k)
                ! higher-order stencil
                ! s(i,hi(2)+1:hi(2)+ng_s,k) = -2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,hi(2)+1:hi(2)+ng_s,k) = -s(i,hi(2),k)
                ! higher-order stencil
                ! s(i,hi(2)+1:hi(2)+ng_s,k) = -2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3)-ng_s,hi(3)+ng_s
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do k=lo(3),hi(3)+1
             do i=lo(1)-ng_s,hi(1)+ng_s
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,lo(3)-ng_s:lo(3)-1) = -s(i,j,lo(3))
                ! higher-order stencil
                ! s(i,j,lo(3)-ng_s:lo(3)-1) = -2.d0*s(i,j,lo(3)) + (1.d0/3.d0)*s(i,j,lo(3)+1)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = -s(i,j,lo(3))
                ! higher-order stencil
                ! s(i,j,lo(3)-ng_s:lo(3)-1) = -2.d0*s(i,j,lo(3)) + (1.d0/3.d0)*s(i,j,lo(3)+1)
             end do
          end do
       else if (bccomp .eq. 3) then 
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(3,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 3) then 
          print *,'physbc_macvel_3d: bc(3,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,hi(3)+1:hi(3)+ng_s) = -s(i,j,hi(3))
                ! higher-order stencil
                ! s(i,j,hi(3)+1:hi(3)+ng_s) = -2.d0*s(i,j,hi(3)) + (1.d0/3.d0)*s(i,j,hi(3)-1)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = -s(i,j,hi(3))
                ! higher-order stencil
                ! s(i,j,hi(3)+1:hi(3)+ng_s) = -2.d0*s(i,j,hi(3)) + (1.d0/3.d0)*s(i,j,hi(3)-1)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel
       end if
    else if (bc(3,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1),hi(1)+1
             do j=lo(2)-ng_s,hi(2)+ng_s
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          ! since du_n/dt=0, the bc is du_t/dn=0
          do i=lo(1)-ng_s,hi(1)+ng_s
             do j=lo(2),hi(2)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3))
             end do
          end do
       else if (bccomp .eq. 3) then
          print *,'physbc_macvel_3d: bc(3,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_3d

  subroutine physbc_macvel_3d_inhomogeneous(s,ng_s,v_bc1,v_bc2,ng_v,lo,hi,bc,bccomp,dx)

    integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_v
    real(kind=dp_t), intent(inout) ::     s(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
    real(kind=dp_t), intent(in   ) :: v_bc1(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t), intent(in   ) :: v_bc2(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    integer        , intent(in   ) :: bc(:,:)
    integer        , intent(in   ) :: bccomp
    real(kind=dp_t), intent(in   ) :: dx(:)

    ! if bccomp=1, v_bc1 represents x-vel on y-faces
    !              v_bc2 represents x-vel on z-faces
    ! if bccomp=2, v_bc1 represents y-vel on x-faces
    !              v_bc2 represents y-vel on z-faces
    ! if bccomp=3, v_bc1 represents z-vel on x-faces
    !              v_bc2 represents z-vel on y-faces

    ! Local variables
    integer :: i,j,k

    if (bccomp .ne. 1 .and. bccomp .ne. 2 .and. bccomp .ne. 3) then
       call bl_error('physbc_macvel_3d_inhomogeneous requires bccomp = 1, 2 or 3')
    end if

!!!!!!!!!!!!!!!!!!
! lo-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = 2.d0*v_bc1(lo(1),j,k) - s(lo(1),j,k)
                ! higher-order stencil
                ! s(lo(1)-ng_s:lo(1)-1,j,k) = (8.d0/3.d0)*v_bc1(lo(1),j,k) - 2.d0*s(lo(1),j,k) + (1.d0/3.d0)*s(lo(1)+1,j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = 2.d0*v_bc1(lo(1),j,k) - s(lo(1),j,k)
                ! higher-order stencil
                ! s(lo(1)-ng_s:lo(1)-1,j,k) = (8.d0/3.d0)*v_bc1(lo(1),j,k) - 2.d0*s(lo(1),j,k) + (1.d0/3.d0)*s(lo(1)+1,j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d_inhomogeneous: bc(1,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k) - dx(1)*v_bc1(lo(1),j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                s(lo(1)-ng_s:lo(1)-1,j,k) = s(lo(1),j,k) - dx(1)*v_bc1(lo(1),j,k)
             end do
          end do
       end if
    else if (bc(1,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(1,1) =',bc(1,1),'for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-x boundary
!!!!!!!!!!!!!!!!!!

    if (bc(1,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = 2.d0*v_bc1(hi(1)+1,j,k) - s(hi(1),j,k)
                ! higher-order stencil
                ! s(hi(1)+1:hi(1)+ng_s,j,k) = (8.d0/3.d0)*v_bc1(hi(1)+1,j,k) - 2.d0*s(hi(1),j,k) + (1.d0/3.d0)*s(hi(1)-1,j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = 2.d0*v_bc1(hi(1)+1,j,k) - s(hi(1),j,k)
                ! higher-order stencil
                ! s(hi(1)+1:hi(1)+ng_s,j,k) = (8.d0/3.d0)*v_bc1(hi(1)+1,j,k) - 2.d0*s(hi(1),j,k) + (1.d0/3.d0)*s(hi(1)-1,j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          print *,'physbc_macvel_3d_inhomogeneous: bc(1,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)+1
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k) + dx(1)*v_bc1(hi(1)+1,j,k)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)+1
             do j=lo(2),hi(2)
                s(hi(1)+1:hi(1)+ng_s,j,k) = s(hi(1),j,k) + dx(1)*v_bc1(hi(1)+1,j,k)
             end do
          end do
       end if
    else if (bc(1,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(1,2) =',bc(1,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = 2.d0*v_bc1(i,lo(2),k) - s(i,lo(2),k)
                ! higher-order stencil
                ! s(i,lo(2)-ng_s:lo(2)-1,k) = (8.d0/3.d0)*v_bc1(i,lo(2),k) - 2.d0*s(i,lo(2),k) + (1.d0/3.d0)*s(i,lo(2)+1,k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1),hi(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = 2.d0*v_bc2(i,lo(2),k) - s(i,lo(2),k)
                ! higher-order stencil
                ! s(i,lo(2)-ng_s:lo(2)-1,k) = (8.d0/3.d0)*v_bc2(i,lo(2),k) - 2.d0*s(i,lo(2),k) + (1.d0/3.d0)*s(i,lo(2)+1,k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)+1
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k) - dx(2)*v_bc1(i,lo(2),k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d_inhomogeneous: bc(2,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1),hi(1)
                s(i,lo(2)-ng_s:lo(2)-1,k) = s(i,lo(2),k) - dx(2)*v_bc2(i,lo(2),k)
             end do
          end do
       end if
    else if (bc(2,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(2,1) =',bc(2,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-y boundary
!!!!!!!!!!!!!!!!!!

    if (bc(2,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = 2.d0*v_bc1(i,hi(2)+1,k) - s(i,hi(2),k)
                ! higher-order stencil
                ! s(i,hi(2)+1:hi(2)+ng_s,k) = (8.d0/3.d0)*v_bc1(i,hi(2)+1,k) - 2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1),hi(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = 2.d0*v_bc2(i,hi(2)+1,k) - s(i,hi(2),k)
                ! higher-order stencil
                ! s(i,hi(2)+1:hi(2)+ng_s,k) = (8.d0/3.d0)*v_bc2(i,hi(2)+1,k) - 2.d0*s(i,hi(2),k) + (1.d0/3.d0)*s(i,hi(2)-1,k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)
             do i=lo(1),hi(1)+1
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k) + dx(2)*v_bc1(i,hi(2)+1,k)
             end do
          end do
       else if (bccomp .eq. 2) then
          print *,'physbc_macvel_3d_inhomogeneous: bc(2,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       else if (bccomp .eq. 3) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do k=lo(3),hi(3)+1
             do i=lo(1),hi(1)
                s(i,hi(2)+1:hi(2)+ng_s,k) = s(i,hi(2),k) + dx(2)*v_bc2(i,hi(2)+1,k)
             end do
          end do
       end if
    else if (bc(2,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(2,2) =',bc(2,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! lo-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,1) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = 2.d0*v_bc2(i,j,lo(3)) - s(i,j,lo(3))
                ! higher-order stencil
                ! s(i,j,lo(3)-ng_s:lo(3)-1) = (8.d0/3.d0)*v_bc2(i,j,lo(3)) - 2.d0*s(i,j,lo(3)) + (1.d0/3.d0)*s(i,j,lo(3)+1)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = 2.d0*v_bc2(i,j,lo(3)) - s(i,j,lo(3))
                ! higher-order stencil
                ! s(i,j,lo(3)-ng_s:lo(3)-1) = (8.d0/3.d0)*v_bc2(i,j,lo(3)) - 2.d0*s(i,j,lo(3)) + (1.d0/3.d0)*s(i,j,lo(3)+1)
             end do
          end do
       else if (bccomp .eq. 3) then 
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(3,1) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3)) - dx(3)*v_bc2(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                s(i,j,lo(3)-ng_s:lo(3)-1) = s(i,j,lo(3)) - dx(3)*v_bc2(i,j,lo(3))
             end do
          end do
       else if (bccomp .eq. 3) then 
          print *,'physbc_macvel_3d_inhomogeneous: bc(3,1) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,1) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(3,1) =',bc(3,1),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

!!!!!!!!!!!!!!!!!!
! hi-z boundary
!!!!!!!!!!!!!!!!!!

    if (bc(3,2) .eq. DIR_VEL) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = 2.d0*v_bc2(i,j,hi(3)+1) - s(i,j,hi(3))
                ! higher-order stencil
                ! s(i,j,hi(3)+1:hi(3)+ng_s) = (8.d0/3.d0)*v_bc2(i,j,hi(3)+1) - 2.d0*s(i,j,hi(3)) + (1.d0/3.d0)*s(i,j,hi(3)-1)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet velocity condition at boundary
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = 2.d0*v_bc2(i,j,hi(3)+1) - s(i,j,hi(3))
                ! higher-order stencil
                ! s(i,j,hi(3)+1:hi(3)+ng_s) = (8.d0/3.d0)*v_bc2(i,j,hi(3)+1) - 2.d0*s(i,j,hi(3)) + (1.d0/3.d0)*s(i,j,hi(3)-1)
             end do
          end do
       else if (bccomp .eq. 3) then
          ! normal velocity
          ! shouldn't have to do anything; this case is covered in physbc_domainvel_inhomogeneous
       end if
    else if (bc(3,2) .eq. DIR_TRACT) then
       if (bccomp .eq. 1) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3)) + dx(3)*v_bc2(i,j,hi(3)+1)
             end do
          end do
       else if (bccomp .eq. 2) then
          ! transverse velocity
          ! two point stencil using homogeneous dirichlet traction condition at boundary
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                s(i,j,hi(3)+1:hi(3)+ng_s) = s(i,j,hi(3)) + dx(3)*v_bc2(i,j,hi(3)+1)
             end do
          end do
       else if (bccomp .eq. 3) then
          print *,'physbc_macvel_3d_inhomogeneous: bc(3,2) = DIR_TRACT for bccomp =',bccomp
          call bl_error('NOT SUPPORTED')
       end if
    else if (bc(3,2) .eq. INTERIOR) then
       ! either periodic or interior; do nothing
    else
       print *,'physbc_macvel_3d_inhomogeneous: bc(3,2) =',bc(3,2),' for bccomp =',bccomp
       call bl_error('NOT SUPPORTED')
    end if

  end subroutine physbc_macvel_3d_inhomogeneous

  subroutine set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time,the_bc_level)

    ! vel_bc_n(nlevs,dm) are the normal velocities
    ! in 2D, vel_bc_t(nlevs,2) respresents
    !   1. y-velocity bc on x-faces (nodal)
    !   2. x-velocity bc on y-faces (nodal)
    ! in 3D, vel_bc_t(nlevs,6) represents
    !   1. y-velocity bc on x-faces (nodal in y and x)
    !   2. z-velocity bc on x-faces (nodal in z and x)
    !   3. x-velocity bc on y-faces (nodal in x and y)
    !   4. z-velocity bc on y-faces (nodal in z and y)
    !   5. x-velocity bc on z-faces (nodal in x and z)
    !   6. y-velocity bc on z-faces (nodal in y and z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)   ! nodal in 2D, edge-based in 3D
    real(kind=dp_t), intent(in   ) :: dx(:,:),time
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,nlevs,i,dm,ng_n,ng_t,ng_e
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: nxp(:,:,:,:)
    real(kind=dp_t), pointer :: nyp(:,:,:,:)
    real(kind=dp_t), pointer :: nzp(:,:,:,:)
    real(kind=dp_t), pointer :: t1p(:,:,:,:)
    real(kind=dp_t), pointer :: t2p(:,:,:,:)
    real(kind=dp_t), pointer :: t3p(:,:,:,:)
    real(kind=dp_t), pointer :: t4p(:,:,:,:)
    real(kind=dp_t), pointer :: t5p(:,:,:,:)
    real(kind=dp_t), pointer :: t6p(:,:,:,:)
    real(kind=dp_t), pointer :: enp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:)
    real(kind=dp_t), pointer :: ep2(:,:,:,:)
    real(kind=dp_t), pointer :: ep3(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"set_inhomogeneous_vel_bcs")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = vel_bc_n(1,1)%ng
    ng_t = vel_bc_t(1,1)%ng
    ng_e = eta_ed(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(vel_bc_n(n,1))
          nxp => dataptr(vel_bc_n(n,1),i)
          nyp => dataptr(vel_bc_n(n,2),i)
          t1p => dataptr(vel_bc_t(n,1),i)
          t2p => dataptr(vel_bc_t(n,2),i)
          enp => dataptr(eta_ed(n,1),i)
          lo = lwb(get_box(vel_bc_n(n,1),i))
          hi = upb(get_box(vel_bc_n(n,1),i))
          select case (dm)
          case (2)
             call set_inhomogeneous_vel_bcs_2d(nxp(:,:,1,1),nyp(:,:,1,1),ng_n, &
                                               t1p(:,:,1,1),t2p(:,:,1,1),ng_t, &
                                               enp(:,:,1,1),ng_e, &
                                               lo,hi,dx(n,:),time, &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             nzp => dataptr(vel_bc_n(n,3),i)
             t3p => dataptr(vel_bc_t(n,3),i)
             t4p => dataptr(vel_bc_t(n,4),i)
             t5p => dataptr(vel_bc_t(n,5),i)
             t6p => dataptr(vel_bc_t(n,6),i)
             ep1 => dataptr(eta_ed(n,1),i)
             ep2 => dataptr(eta_ed(n,2),i)
             ep3 => dataptr(eta_ed(n,3),i)
             call set_inhomogeneous_vel_bcs_3d(nxp(:,:,:,1),nyp(:,:,:,1),nzp(:,:,:,1),ng_n, &
                                               t1p(:,:,:,1),t2p(:,:,:,1),t3p(:,:,:,1), &
                                               t4p(:,:,:,1),t5p(:,:,:,1),t6p(:,:,:,1),ng_t, &
                                               ep1(:,:,:,1),ep2(:,:,:,1),ep3(:,:,:,1),ng_e, &
                                               lo,hi,dx(n,:),time, &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          end select
       end do
    end do
    
    call destroy(bpt)

  end subroutine set_inhomogeneous_vel_bcs

  subroutine set_inhomogeneous_vel_bcs_2d(vel_bc_nx,vel_bc_ny,ng_n, &
                                          vel_bc_tyx,vel_bc_txy,ng_t, &
                                          eta_nd,ng_e,lo,hi,dx,time,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t,ng_e
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) ::     eta_nd(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: x,y

    !!!!!!!!!!!!!!!!!!!!
    ! normal velocities
    !!!!!!!!!!!!!!!!!!!!

    ! xvel, lo x-faces
    if (bc(1,1,1) .eq. DIR_VEL) then
       x = prob_lo(1)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          vel_bc_nx(lo(1),j) = inhomogeneous_bc_val_2d(1,x,y,time)
       end do
    end if

    ! xvel, hi x-faces
    if (bc(1,2,1) .eq. DIR_VEL) then
       x = prob_hi(1)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          vel_bc_nx(hi(1)+1,j) = inhomogeneous_bc_val_2d(1,x,y,time)
       end do
    end if

    ! yvel, lo y-faces
    if (bc(2,1,2) .eq. DIR_VEL) then
       y = prob_lo(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          vel_bc_ny(i,lo(2)) = inhomogeneous_bc_val_2d(2,x,y,time)
       end do
    end if

    ! yvel, hi y-faces
    if (bc(2,2,2) .eq. DIR_VEL) then
       y = prob_hi(2)
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
          vel_bc_ny(i,hi(2)+1) = inhomogeneous_bc_val_2d(2,x,y,time)
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! yvel, lo x-faces
    if (bc(1,1,2) .eq. DIR_VEL .or. bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(lo(1),j) = inhomogeneous_bc_val_2d(2,x,y,time)
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,1,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          vel_bc_tyx(lo(1),j) = vel_bc_tyx(lo(1),j) / eta_nd(lo(1),j)
       end do
    end if

    ! yvel, hi x-faces
    if (bc(1,2,2) .eq. DIR_VEL .or. bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          vel_bc_tyx(hi(1)+1,j) = inhomogeneous_bc_val_2d(2,x,y,time)
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,2,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          vel_bc_tyx(hi(1)+1,j) = vel_bc_tyx(hi(1)+1,j) / eta_nd(hi(1)+1,j)
       end do
    end if

    ! xvel, lo y-faces
    if (bc(2,1,1) .eq. DIR_VEL .or. bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,lo(2)) = inhomogeneous_bc_val_2d(1,x,y,time)
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,1,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          vel_bc_txy(i,lo(2)) = vel_bc_txy(i,lo(2)) / eta_nd(i,lo(2))
       end do
    end if

    ! xvel, hi y-faces
    if (bc(2,2,1) .eq. DIR_VEL .or. bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          vel_bc_txy(i,hi(2)+1) = inhomogeneous_bc_val_2d(1,x,y,time)
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,2,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          vel_bc_txy(i,hi(2)+1) = vel_bc_txy(i,hi(2)+1) / eta_nd(i,hi(2)+1)
       end do
    end if

  end subroutine set_inhomogeneous_vel_bcs_2d

  subroutine set_inhomogeneous_vel_bcs_3d(vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_n, &
                                          vel_bc_tyx,vel_bc_tzx,vel_bc_txy,vel_bc_tzy, &
                                          vel_bc_txz,vel_bc_tyz,ng_t, &
                                          eta_xy,eta_xz,eta_yz,ng_e,lo,hi,dx,time,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t,ng_e
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_nz(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(in   ) ::     eta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::     eta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::     eta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: dx(:),time
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: x,y,z

    !!!!!!!!!!!!!!!!!!!!
    ! normal velocities
    !!!!!!!!!!!!!!!!!!!!

    ! xvel, lo x-faces
    if (bc(1,1,1) .eq. DIR_VEL) then
       x = prob_lo(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_nx(lo(1),j,k) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if

    ! xvel, hi x-faces
    if (bc(1,2,1) .eq. DIR_VEL) then
       x = prob_hi(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_nx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if
    
    ! yvel, lo y-faces
    if (bc(2,1,2) .eq. DIR_VEL) then
       y = prob_lo(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_ny(i,lo(2),k) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if

    ! yvel, hi y-faces
    if (bc(2,2,2) .eq. DIR_VEL) then
       y = prob_hi(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_ny(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if
    
    ! zvel, lo z-faces
    if (bc(3,1,3) .eq. DIR_VEL) then
       z = prob_lo(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_nz(i,j,lo(3)) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    ! zvel, hi z-faces
    if (bc(3,2,3) .eq. DIR_VEL) then
       z = prob_hi(3)
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_nz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!
    
    ! yvel, lo x-faces
    if (bc(1,1,2) .eq. DIR_VEL .or. bc(1,1,2) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(lo(1),j,k) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,1,2) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             vel_bc_tyx(lo(1),j,k) = vel_bc_tyx(lo(1),j,k) / eta_xy(lo(1),j,k)
          end do
       end do
    end if

    ! yvel, hi x-faces
    if (bc(1,2,2) .eq. DIR_VEL .or. bc(1,2,2) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do j=lo(2),hi(2)+1
             y = prob_lo(2) + dble(j)*dx(2)
             vel_bc_tyx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,2,2) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             vel_bc_tyx(hi(1)+1,j,k) = vel_bc_tyx(hi(1)+1,j,k) / eta_xy(hi(1)+1,j,k)
          end do
       end do
    end if

    ! zvel, lo x-faces
    if (bc(1,1,3) .eq. DIR_VEL .or. bc(1,1,3) .eq. DIR_TRACT) then
       x = prob_lo(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(lo(1),j,k) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,1,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             vel_bc_tzx(lo(1),j,k) = vel_bc_tzx(lo(1),j,k) / eta_xz(lo(1),j,k)
          end do
       end do
    end if

    ! zvel, hi x-faces
    if (bc(1,2,3) .eq. DIR_VEL .or. bc(1,2,3) .eq. DIR_TRACT) then
       x = prob_hi(1)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_tzx(hi(1)+1,j,k) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(1,2,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             vel_bc_tzx(hi(1)+1,j,k) = vel_bc_tzx(hi(1)+1,j,k) / eta_xz(hi(1)+1,j,k)
          end do
       end do
    end if

    ! xvel, lo y-faces
    if (bc(2,1,1) .eq. DIR_VEL .or. bc(2,1,1) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,lo(2),k) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,1,1) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)+1
             vel_bc_txy(i,lo(2),k) = vel_bc_txy(i,lo(2),k) / eta_xy(i,lo(2),k)
          end do
       end do
    end if

    ! xvel, hi y-faces
    if (bc(2,2,1) .eq. DIR_VEL .or. bc(2,2,1) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+0.5d0)*dx(3)
          do i=lo(1),hi(1)+1
             x = prob_lo(1) + dble(i)*dx(1)
             vel_bc_txy(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,2,1) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)+1
             vel_bc_txy(i,hi(2)+1,k) = vel_bc_txy(i,hi(2)+1,k) / eta_xy(i,hi(2)+1,k)
          end do
       end do
    end if

    ! zvel, lo y-faces
    if (bc(2,1,3) .eq. DIR_VEL .or. bc(2,1,3) .eq. DIR_TRACT) then
       y = prob_lo(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,lo(2),k) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,1,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do i=lo(1),hi(1)
             vel_bc_tzy(i,lo(2),k) = vel_bc_tzy(i,lo(2),k) / eta_yz(i,lo(2),k)
          end do
       end do
    end if

    ! zvel, hi y-faces
    if (bc(2,2,3) .eq. DIR_VEL .or. bc(2,2,3) .eq. DIR_TRACT) then
       y = prob_hi(2)
       do k=lo(3),hi(3)+1
          z = prob_lo(3) + dble(k)*dx(3)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tzy(i,hi(2)+1,k) = inhomogeneous_bc_val_3d(3,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(2,2,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do i=lo(1),hi(1)
             vel_bc_tzy(i,hi(2)+1,k) = vel_bc_tzy(i,hi(2)+1,k) / eta_yz(i,hi(2)+1,k)
          end do
       end do
    end if

    ! xvel, lo z-faces
    if (bc(3,1,1) .eq. DIR_VEL .or. bc(3,1,1) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,lo(3)) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(3,1,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          do j=lo(2),hi(2)
             vel_bc_txz(i,j,lo(3)) = vel_bc_txz(i,j,lo(3)) / eta_xz(i,j,lo(3))
          end do
       end do
    end if

    ! xvel, hi z-faces
    if (bc(3,2,1) .eq. DIR_VEL .or. bc(3,2,1) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do i=lo(1),hi(1)+1
          x = prob_lo(1) + dble(i)*dx(1)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+0.5d0)*dx(2)
             vel_bc_txz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(1,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(3,2,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          do j=lo(2),hi(2)
             vel_bc_txz(i,j,hi(3)+1) = vel_bc_txz(i,j,hi(3)+1) / eta_xz(i,j,hi(3)+1)
          end do
       end do
    end if

    ! yvel, lo z-faces
    if (bc(3,1,2) .eq. DIR_VEL .or. bc(3,1,2) .eq. DIR_TRACT) then
       z = prob_lo(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,lo(3)) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(3,1,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             vel_bc_tyz(i,j,lo(3)) = vel_bc_tyz(i,j,lo(3)) / eta_yz(i,j,lo(3))
          end do
       end do
    end if

    ! yvel, hi z-faces
    if (bc(3,2,2) .eq. DIR_VEL .or. bc(3,2,2) .eq. DIR_TRACT) then
       z = prob_hi(3)
       do j=lo(2),hi(2)+1
          y = prob_lo(2) + dble(j)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
             vel_bc_tyz(i,j,hi(3)+1) = inhomogeneous_bc_val_3d(2,x,y,z,time)
          end do
       end do
    end if

    ! divide by eta for traction condition
    if (bc(3,2,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             vel_bc_tyz(i,j,hi(3)+1) = vel_bc_tyz(i,j,hi(3)+1) / eta_yz(i,j,hi(3)+1)
          end do
       end do
    end if

  end subroutine set_inhomogeneous_vel_bcs_3d

  subroutine modify_traction_bcs(mla,vel_bc_n,vel_bc_t,dx,the_bc_level)

    ! This subroutine modifies the traction bc, f / eta, by subtracting
    ! off the dv_n/dt part for v_t

    ! vel_bc_n(nlevs,dm) are the normal velocities
    ! in 2D, vel_bc_t(nlevs,2) respresents
    !   1. y-velocity bc on x-faces (nodal)
    !   2. x-velocity bc on y-faces (nodal)
    ! in 3D, vel_bc_t(nlevs,6) represents
    !   1. y-velocity bc on x-faces (nodal in y and x)
    !   2. z-velocity bc on x-faces (nodal in z and x)
    !   3. x-velocity bc on y-faces (nodal in x and y)
    !   4. z-velocity bc on y-faces (nodal in z and y)
    !   5. x-velocity bc on z-faces (nodal in x and z)
    !   6. y-velocity bc on z-faces (nodal in y and z)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: vel_bc_n(:,:)
    type(multifab) , intent(inout) :: vel_bc_t(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: n,nlevs,i,dm,ng_n,ng_t
    integer :: lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: nxp(:,:,:,:)
    real(kind=dp_t), pointer :: nyp(:,:,:,:)
    real(kind=dp_t), pointer :: nzp(:,:,:,:)
    real(kind=dp_t), pointer :: t1p(:,:,:,:)
    real(kind=dp_t), pointer :: t2p(:,:,:,:)
    real(kind=dp_t), pointer :: t3p(:,:,:,:)
    real(kind=dp_t), pointer :: t4p(:,:,:,:)
    real(kind=dp_t), pointer :: t5p(:,:,:,:)
    real(kind=dp_t), pointer :: t6p(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"modify_traction_bcs")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_n = vel_bc_n(1,1)%ng
    ng_t = vel_bc_t(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(vel_bc_n(n,1))
          nxp => dataptr(vel_bc_n(n,1),i)
          nyp => dataptr(vel_bc_n(n,2),i)
          t1p => dataptr(vel_bc_t(n,1),i)
          t2p => dataptr(vel_bc_t(n,2),i)
          lo = lwb(get_box(vel_bc_n(n,1),i))
          hi = upb(get_box(vel_bc_n(n,1),i))
          select case (dm)
          case (2)
             call modify_traction_bcs_2d(nxp(:,:,1,1),nyp(:,:,1,1),ng_n, &
                                               t1p(:,:,1,1),t2p(:,:,1,1),ng_t, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))
          case (3)
             nzp => dataptr(vel_bc_n(n,3),i)
             t3p => dataptr(vel_bc_t(n,3),i)
             t4p => dataptr(vel_bc_t(n,4),i)
             t5p => dataptr(vel_bc_t(n,5),i)
             t6p => dataptr(vel_bc_t(n,6),i)
             call modify_traction_bcs_3d(nxp(:,:,:,1),nyp(:,:,:,1),nzp(:,:,:,1),ng_n, &
                                               t1p(:,:,:,1),t2p(:,:,:,1),t3p(:,:,:,1), &
                                               t4p(:,:,:,1),t5p(:,:,:,1),t6p(:,:,:,1),ng_t, &
                                               lo,hi,dx(n,:), &
                                               the_bc_level(n)%adv_bc_level_array(i,:,:,:))

          end select
       end do
    end do
    
    call destroy(bpt)

  end subroutine modify_traction_bcs

  subroutine modify_traction_bcs_2d(vel_bc_nx,vel_bc_ny,ng_n, &
                                          vel_bc_tyx,vel_bc_txy,ng_t,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j
    real(kind=dp_t) :: dxinv

    dxinv = 1.d0/dx(1)

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!

    ! yvel, lo x-faces
    ! subtract dvx/dy
    if (bc(1,1,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          vel_bc_tyx(lo(1),j) = vel_bc_tyx(lo(1),j) &
               - (vel_bc_nx(lo(1),j)-vel_bc_nx(lo(1),j-1)) * dxinv
       end do
    end if

    ! yvel, hi x-faces
    ! subtract dvx/dy
    if (bc(1,2,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          vel_bc_tyx(hi(1)+1,j) = vel_bc_tyx(hi(1)+1,j) &
               - (vel_bc_nx(hi(1)+1,j)-vel_bc_nx(hi(1)+1,j-1)) * dxinv
       end do
    end if

    ! xvel, lo y-faces
    ! subtract dvy/dx
    if (bc(2,1,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          vel_bc_txy(i,lo(2)) = vel_bc_txy(i,lo(2)) &
               - (vel_bc_ny(i,lo(2))-vel_bc_ny(i-1,lo(2))) * dxinv

       end do
    end if

    ! xvel, hi y-faces
    ! subtract dvy/dx
    if (bc(2,2,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          vel_bc_txy(i,hi(2)+1) = vel_bc_txy(i,hi(2)+1) &
               - (vel_bc_ny(i,hi(2)+1)-vel_bc_ny(i-1,hi(2)+1)) * dxinv
       end do
    end if

  end subroutine modify_traction_bcs_2d

  subroutine modify_traction_bcs_3d(vel_bc_nx,vel_bc_ny,vel_bc_nz,ng_n, &
                                          vel_bc_tyx,vel_bc_tzx,vel_bc_txy,vel_bc_tzy, &
                                          vel_bc_txz,vel_bc_tyz,ng_t,lo,hi,dx,bc)

    integer        , intent(in   ) :: lo(:),hi(:),ng_n,ng_t
    real(kind=dp_t), intent(inout) ::  vel_bc_nx(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_ny(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) ::  vel_bc_nz(lo(1)-ng_n:,lo(2)-ng_n:,lo(3)-ng_n:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzx(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tzy(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_txz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(inout) :: vel_bc_tyz(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:,:)

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: dxinv

    dxinv = 1.d0/dx(1)

    !!!!!!!!!!!!!!!!!!!!!!!!
    ! transverse velocities
    !!!!!!!!!!!!!!!!!!!!!!!!
    
    ! yvel, lo x-faces
    ! subtract dvx/dy
    if (bc(1,1,2) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             vel_bc_tyx(lo(1),j,k) = vel_bc_tyx(lo(1),j,k) &
                  - (vel_bc_nx(lo(1),j,k)-vel_bc_nx(lo(1),j-1,k)) * dxinv
          end do
       end do

    end if

    ! yvel, hi x-faces
    ! subtract dvx/dy
    if (bc(1,2,2) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             vel_bc_tyx(hi(1)+1,j,k) = vel_bc_tyx(hi(1)+1,j,k) &
                  - (vel_bc_nx(hi(1)+1,j,k)-vel_bc_nx(hi(1)+1,j-1,k)) * dxinv
          end do
       end do

    end if

    ! zvel, lo x-faces
    ! subtract dvx/dz
    if (bc(1,1,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             vel_bc_tzx(lo(1),j,k) = vel_bc_tzx(lo(1),j,k) &
                  - (vel_bc_nx(lo(1),j,k)-vel_bc_nx(lo(1),j,k-1)) * dxinv
          end do
       end do

    end if

    ! zvel, hi x-faces
    ! subtract dvx/dz
    if (bc(1,2,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             vel_bc_tzx(hi(1)+1,j,k) = vel_bc_tzx(hi(1)+1,j,k) &
                  - (vel_bc_nx(hi(1)+1,j,k)-vel_bc_nx(hi(1)+1,j,k-1)) * dxinv
          end do
       end do

    end if

    ! xvel, lo y-faces
    ! subtract dvy/dx
    if (bc(2,1,1) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)+1
             vel_bc_txy(i,lo(2),k) = vel_bc_txy(i,lo(2),k) &
                  - (vel_bc_ny(i,lo(2),k)-vel_bc_ny(i-1,lo(2),k)) * dxinv
          end do
       end do

    end if

    ! xvel, hi y-faces
    ! subtract dvy/dx
    if (bc(2,2,1) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)
          do i=lo(1),hi(1)+1
             vel_bc_txy(i,hi(2)+1,k) = vel_bc_txy(i,hi(2)+1,k) &
                  - (vel_bc_ny(i,hi(2)+1,k)-vel_bc_ny(i-1,hi(2)+1,k)) * dxinv
          end do
       end do

    end if

    ! zvel, lo y-faces
    ! subtract dvy/dz
    if (bc(2,1,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do i=lo(1),hi(1)
             vel_bc_tzy(i,lo(2),k) = vel_bc_tzy(i,lo(2),k) &
                  - (vel_bc_ny(i,lo(2),k)-vel_bc_ny(i,lo(2),k-1)) * dxinv
          end do
       end do

    end if

    ! zvel, hi y-faces
    ! subtract dvy/dz
    if (bc(2,2,3) .eq. DIR_TRACT) then
       do k=lo(3),hi(3)+1
          do i=lo(1),hi(1)
             vel_bc_tzy(i,hi(2)+1,k) = vel_bc_tzy(i,hi(2)+1,k) &
                  - (vel_bc_ny(i,hi(2)+1,k)-vel_bc_ny(i,hi(2)+1,k-1)) * dxinv
          end do
       end do

    end if

    ! xvel, lo z-faces
    ! subtract dvz/dx
    if (bc(3,1,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          do j=lo(2),hi(2)
             vel_bc_txz(i,j,lo(3)) = vel_bc_txz(i,j,lo(3)) &
                  - (vel_bc_nz(i,j,lo(3))-vel_bc_nz(i-1,j,lo(3))) * dxinv
          end do
       end do

    end if

    ! xvel, hi z-faces
    ! subtract dvz/dx
    if (bc(3,2,1) .eq. DIR_TRACT) then
       do i=lo(1),hi(1)+1
          do j=lo(2),hi(2)
             vel_bc_txz(i,j,hi(3)+1) = vel_bc_txz(i,j,hi(3)+1) &
                  - (vel_bc_nz(i,j,hi(3)+1)-vel_bc_nz(i-1,j,hi(3)+1)) * dxinv
          end do
       end do

    end if

    ! yvel, lo z-faces
    ! subtract dvz/dy
    if (bc(3,1,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             vel_bc_tyz(i,j,lo(3)) = vel_bc_tyz(i,j,lo(3)) &
                  - (vel_bc_nz(i,j,lo(3))-vel_bc_nz(i,j-1,lo(3))) * dxinv
          end do
       end do

    end if

    ! yvel, hi z-faces
    ! subtract dvz/dy
    if (bc(3,2,2) .eq. DIR_TRACT) then
       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             vel_bc_tyz(i,j,hi(3)+1) = vel_bc_tyz(i,j,hi(3)+1) &
                  - (vel_bc_nz(i,j,hi(3)+1)-vel_bc_nz(i,j-1,hi(3)+1)) * dxinv
          end do
       end do

    end if

  end subroutine modify_traction_bcs_3d

end module multifab_physbc_stag_module
