module estdt_module 

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: cfl

  implicit none

  private

  public :: estdt

contains

  subroutine estdt (mla ,umac, dx, dt)
    
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(  out) :: dt

    real(kind=dp_t), pointer:: ump(:,:,:,:), vmp(:,:,:,:), wmp(:,:,:,:)

    integer lo(mla%dim), hi(mla%dim), ng, i, n, nlevs, dm
    real(kind=dp_t) :: dt_proc, dt_grid

    type(bl_prof_timer),save :: bpt

    call build(bpt,"estdt")

    ng = umac(1,1)%ng
    dm = mla%dim
    nlevs = mla%nlevel

    dt_proc  = 1.d20

    do n=1,nlevs
       do i=1,nfabs(umac(n,1))
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          lo =  lwb(get_box(umac(n,1), i))
          hi =  upb(get_box(umac(n,1), i))
          dt_grid = 1.d20
          select case (dm)
          case (2)
             call estdt_2d(ump(:,:,1,1), vmp(:,:,1,1), ng, &
                           lo, hi, dx(n,:), dt_grid)
          case (3)
             wmp => dataptr(umac(n,3), i)
             call estdt_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng, &
                           lo, hi, dx(n,:), dt_grid)
          end select
          dt_proc = min(dt_grid, dt_proc)
       end do
    end do

    ! This sets dt to be the min of dt_proc over all processors.
    call parallel_reduce(dt, dt_proc, MPI_MIN)

    dt = dt * cfl

    call destroy(bpt)

  end subroutine estdt

  subroutine estdt_2d(umac,vmac,ng,lo,hi,dx,dt)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: umac(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: vmac(lo(1)-ng:,lo(2)-ng:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    ! Local variables
    real (kind = dp_t)  u,v,eps
    integer :: i, j

    u  = 0.d0
    v  = 0.d0

    eps = 1.d-10

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
       u = max(u,abs(umac(i,j)))
    enddo
    enddo

    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
       v = max(v,abs(vmac(i,j)))
    enddo
    enddo

    if (u .gt. eps) dt = min(dt,dx(1)/u)
    if (v .gt. eps) dt = min(dt,dx(2)/v)

  end subroutine estdt_2d

  subroutine estdt_3d(umac,vmac,wmac,ng,lo,hi,dx,dt)

    integer           , intent(in   ) :: lo(:), hi(:), ng
    real (kind = dp_t), intent(in   ) :: umac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: vmac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: wmac(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)
    real (kind = dp_t), intent(in   ) :: dx(:)
    real (kind = dp_t), intent(inout) :: dt

    ! Local variables
    real (kind = dp_t)  u,v,w,eps
    integer :: i, j, k

    u  = 0.d0
    v  = 0.d0
    w = 0.d0

    eps = 1.d-10

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
       u = max(u,abs(umac(i,j,k)))
    enddo
    enddo
    enddo

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
       v = max(v,abs(vmac(i,j,k)))
    enddo
    enddo
    enddo

    do k=lo(3),hi(3)+1
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       w = max(w,abs(wmac(i,j,k)))
    enddo
    enddo
    enddo

    if (u .gt. eps) dt = min(dt,dx(1)/u)
    if (v .gt. eps) dt = min(dt,dx(2)/v)
    if (w .gt. eps) dt = min(dt,dx(3)/w)

  end subroutine estdt_3d

end module estdt_module
