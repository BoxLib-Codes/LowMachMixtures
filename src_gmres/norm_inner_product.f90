 module norm_inner_product_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use probin_common_module, only: total_volume

  implicit none

  private
  public :: stag_inner_prod, cc_inner_prod, stag_l2_norm, cc_l2_norm, &
       stag_l1_norm, sum_umac_press

contains

  subroutine stag_inner_prod(mla,m1,comp1,m2,comp2,prod_val)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m1(:,:), m2(:,:)  ! (mla%nlevel,mla%dim)
    integer        , intent(in   ) :: comp1, comp2
    real(kind=dp_t), intent(inout) :: prod_val(:)

    ! local
    integer :: i,comp,dm,nlevs,ng_m1,ng_m2
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: inner_prod_proc(mla%dim), inner_prod_grid(mla%dim)

    real(kind=dp_t), pointer :: m1xp(:,:,:,:),m1yp(:,:,:,:),m1zp(:,:,:,:)
    real(kind=dp_t), pointer :: m2xp(:,:,:,:),m2yp(:,:,:,:),m2zp(:,:,:,:)
    
    type(multifab) :: temp_cc
    
    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stag_inner_prod")

    call multifab_build(temp_cc,mla%la(1),1,0)

    inner_prod_proc = 0.d0
    nlevs = mla%nlevel
    dm = mla%dim

    ng_m1 = m1(1,1)%ng
    ng_m2 = m2(1,1)%ng

    if (nlevs .gt. 1) then
       call bl_error('stag_inner_prod not written for multilevel yet')
    end if

    !$omp parallel private(mfi,i,xnodalbox,xlo,xhi,ynodalbox,ylo,yhi) &
    !$omp private(znodalbox,zlo,zhi,m1xp,m2xp,m1yp,m2yp,m1zp,m2zp) &
    !$omp private(lo,hi,inner_prod_grid) &
    !$omp reduction(+:inner_prod_proc)

    call mfiter_build(mfi, temp_cc, tiling=.true.)

      do while (more_tile(mfi))
         i = get_fab_index(mfi)

         xnodalbox = get_nodaltilebox(mfi,1)
         xlo = lwb(xnodalbox)
         xhi = upb(xnodalbox)
         ynodalbox = get_nodaltilebox(mfi,2)
         ylo = lwb(ynodalbox)
         yhi = upb(ynodalbox)
         znodalbox = get_nodaltilebox(mfi,3)
         zlo = lwb(znodalbox)
         zhi = upb(znodalbox)

!    do i=1,nfabs(m1(1,1))
       m1xp => dataptr(m1(1,1), i)   ! for u
       m2xp => dataptr(m2(1,1), i)   ! for v

       m1yp => dataptr(m1(1,2), i)   ! for u
       m2yp => dataptr(m2(1,2), i)   ! for v

       lo = lwb(get_box(m1(1,1), i))
       hi = upb(get_box(m1(1,1), i))
       inner_prod_grid = 0.d0
       select case (dm)
       case (2)
          call stag_inner_prod_2d(m1xp(:,:,1,comp1),m1yp(:,:,1,comp1),ng_m1, &
                                  m2xp(:,:,1,comp2),m2yp(:,:,1,comp2),ng_m2, &
                                  lo, hi, inner_prod_grid,xlo,xhi,ylo,yhi)
       case (3)
          m1zp => dataptr(m1(1,3), i)   ! for w
          m2zp => dataptr(m2(1,3), i)   ! for w
          call stag_inner_prod_3d(m1xp(:,:,:,comp1),m1yp(:,:,:,comp1),m1zp(:,:,:,comp1),ng_m1, &
                                  m2xp(:,:,:,comp2),m2yp(:,:,:,comp2),m2zp(:,:,:,comp2),ng_m2, &
                                  lo, hi, inner_prod_grid,xlo,xhi,ylo,yhi,zlo,zhi)
       end select

       inner_prod_proc(1:dm) = inner_prod_proc(1:dm) + inner_prod_grid(1:dm)

    end do
    !$omp end parallel

!    do comp=1,dm
       call parallel_reduce(prod_val, inner_prod_proc, MPI_SUM)
!    end do

    call multifab_destroy(temp_cc)

    call destroy(bpt)

  end subroutine stag_inner_prod

  subroutine stag_inner_prod_2d(m1x,m1y,ng_m1,m2x,m2y,ng_m2,glo,ghi,inner_prod,xlo,xhi,ylo,yhi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m1, ng_m2
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
    real(kind=dp_t), intent(in   ) :: m1x(glo(1)-ng_m1:,glo(2)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m1y(glo(1)-ng_m1:,glo(2)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m2x(glo(1)-ng_m2:,glo(2)-ng_m2:)
    real(kind=dp_t), intent(in   ) :: m2y(glo(1)-ng_m2:,glo(2)-ng_m2:)
    real(kind=dp_t), intent(inout) :: inner_prod(:)

    ! local
    integer :: i,j

    ! mx, interior cells
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)
             inner_prod(1) = inner_prod(1) + m1x(i,j)*m2x(i,j)
          end do
       end do

    ! mx, boundary cells
       if (xlo(1) .eq. glo(1)) then
       do j=xlo(2),xhi(2)
          inner_prod(1) = inner_prod(1) + 0.5d0*m1x(glo(1),j)*m2x(glo(1),j)
       end do
       end if

       if (xhi(1) .eq. ghi(1)+1) then
       do j=xlo(2),xhi(2)
          inner_prod(1) = inner_prod(1) + 0.5d0*m1x(ghi(1)+1,j)*m2x(ghi(1)+1,j)
       end do
       end if

    ! my, interior cells
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)
             inner_prod(2) = inner_prod(2) + m1y(i,j)*m2y(i,j)
          end do
       end do

    ! my, boundary cells
       if (ylo(2) .eq. glo(2)) then
       do i=ylo(1),yhi(1)
          inner_prod(2) = inner_prod(2) + 0.5d0*m1y(i,glo(2))* m2y(i,glo(2))
       end do
       end if

       if (yhi(2) .eq. ghi(2)+1) then
       do i=ylo(1),yhi(1)
          inner_prod(2) = inner_prod(2) + 0.5d0*m1y(i,ghi(2)+1)*m2y(i,ghi(2)+1)
       end do
       end if

  end subroutine stag_inner_prod_2d

  subroutine stag_inner_prod_3d(m1x,m1y,m1z,ng_m1,m2x,m2y,m2z,ng_m2,glo,ghi,inner_prod, &
                                       xlo,xhi,ylo,yhi,zlo,zhi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m1, ng_m2
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t), intent(in   ) :: m1x(glo(1)-ng_m1:,glo(2)-ng_m1:,glo(3)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m1y(glo(1)-ng_m1:,glo(2)-ng_m1:,glo(3)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m1z(glo(1)-ng_m1:,glo(2)-ng_m1:,glo(3)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m2x(glo(1)-ng_m2:,glo(2)-ng_m2:,glo(3)-ng_m2:)
    real(kind=dp_t), intent(in   ) :: m2y(glo(1)-ng_m2:,glo(2)-ng_m2:,glo(3)-ng_m2:)
    real(kind=dp_t), intent(in   ) :: m2z(glo(1)-ng_m2:,glo(2)-ng_m2:,glo(3)-ng_m2:)
    real(kind=dp_t), intent(inout) :: inner_prod(:)

    ! local
    integer :: i,j,k

    ! x-comp, interior cells
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)
             inner_prod(1) = inner_prod(1)+m1x(i,j,k)*m2x(i,j,k)
          end do
       end do
    end do

    ! x-comp, boundary cells
    if (xlo(1) .eq. glo(1)) then
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          inner_prod(1) = inner_prod(1) + 0.5d0*m1x(glo(1),j,k)*m2x(glo(1),j,k)
       end do
    end do
    end if

    if (xhi(1) .eq. ghi(1)+1) then
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          inner_prod(1) = inner_prod(1) + 0.5d0*m1x(ghi(1)+1,j,k)*m2x(ghi(1)+1,j,k)
       end do
    end do
    end if

    ! y-comp, interior cells
    do k=ylo(3),yhi(3)
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)
             inner_prod(2) = inner_prod(2) + m1y(i,j,k)*m2y(i,j,k)
          end do
       end do
    end do

    ! y-comp, boundary cells
    if (ylo(2) .eq. glo(2)) then
    do k=ylo(3),yhi(3)
       do i=ylo(1),yhi(1)
          inner_prod(2) = inner_prod(2) + 0.5d0*m1y(i,glo(2),k)*m2y(i,glo(2),k)
       end do
    end do
    end if

    if (yhi(2) .eq. ghi(2)+1) then
    do k=ylo(3),yhi(3)
       do i=ylo(1),yhi(1)
          inner_prod(2) = inner_prod(2) + 0.5d0*m1y(i,ghi(2)+1,k)*m2y(i,ghi(2)+1,k)
       end do
    end do
    end if

    ! z-comp, interior cells
    do k=zlo(3),zhi(3)
       do j=zlo(2),zhi(2)
          do i=zlo(1),zhi(1)
             inner_prod(3) = inner_prod(3) + m1z(i,j,k)*m2z(i,j,k)
          end do
       end do
    end do

    ! z-comp, boundary cells
    if (zlo(3) .eq. glo(3)) then
    do j=zlo(2),zhi(2)
       do i=zlo(1),zhi(1)
          inner_prod(3) = inner_prod(3) + 0.5d0*m1z(i,j,glo(3))*m2z(i,j,glo(3))
       end do
    end do
    end if

    if (zhi(3) .eq. ghi(3)+1) then
    do j=zlo(2),zhi(2)
       do i=zlo(1),zhi(1)
          inner_prod(3) = inner_prod(3) + 0.5d0*m1z(i,j,ghi(3)+1)*m2z(i,j,ghi(3)+1)
       end do
    end do
    end if

  end subroutine stag_inner_prod_3d

  subroutine cc_inner_prod(mla,m1,comp1,m2,comp2,prod_val)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m1(:), m2(:)  ! (mla%nlevel)
    integer        , intent(in   ) :: comp1, comp2
    real(kind=dp_t), intent(inout) :: prod_val

    ! local
    integer :: i,dm,nlevs,ng_m1,ng_m2
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: inner_prod_proc, inner_prod_grid

    real(kind=dp_t), pointer :: m1p(:,:,:,:)
    real(kind=dp_t), pointer :: m2p(:,:,:,:)
    
    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"cc_inner_prod")

    inner_prod_proc = 0.d0
    nlevs = mla%nlevel
    dm = mla%dim

    ng_m1 = m1(1)%ng
    ng_m2 = m2(1)%ng

    if (nlevs .gt. 1) then
       call bl_error('cc_inner_prod not written for multilevel yet')
    end if

    !$omp parallel private(mfi,i,tilebox,tlo,thi,m1p,m2p,lo,hi,inner_prod_grid) &
    !$omp reduction(+:inner_prod_proc)

    call mfiter_build(mfi, m1(1), tiling=.true.)

    do while (more_tile(mfi))
       i = get_fab_index(mfi)

       tilebox = get_tilebox(mfi)
       tlo = lwb(tilebox)
       thi = upb(tilebox)
!    do i=1,nfabs(m1(1))
       m1p => dataptr(m1(1), i)   ! for u
       m2p => dataptr(m2(1), i)   ! for v

       lo = lwb(get_box(m1(1), i))
       hi = upb(get_box(m1(1), i))
       inner_prod_grid = 0.d0
       select case (dm)
       case (2)
          call cc_inner_prod_2d(m1p(:,:,1,comp1),ng_m1,m2p(:,:,1,comp2),ng_m2, &
                                lo, hi, inner_prod_grid,tlo,thi)
       case (3)
          call cc_inner_prod_3d(m1p(:,:,:,comp1),ng_m1,m2p(:,:,:,comp2),ng_m2, &
                                lo, hi, inner_prod_grid,tlo,thi)
       end select

       inner_prod_proc = inner_prod_proc + inner_prod_grid

    end do
    !$omp end parallel

    call parallel_reduce(prod_val, inner_prod_proc, MPI_SUM)

    call destroy(bpt)

  end subroutine cc_inner_prod

  subroutine cc_inner_prod_2d(m1,ng_m1,m2,ng_m2,glo,ghi,inner_prod,tlo,thi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m1, ng_m2,tlo(:),thi(:)
    real(kind=dp_t), intent(in   ) :: m1(glo(1)-ng_m1:,glo(2)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m2(glo(1)-ng_m2:,glo(2)-ng_m2:)
    real(kind=dp_t), intent(inout) :: inner_prod

    ! local
    integer :: i,j

    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
          inner_prod = inner_prod + m1(i,j)*m2(i,j)
       end do
    end do

  end subroutine cc_inner_prod_2d

  subroutine cc_inner_prod_3d(m1,ng_m1,m2,ng_m2,glo,ghi,inner_prod,tlo,thi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m1, ng_m2,tlo(:),thi(:)
    real(kind=dp_t), intent(in   ) :: m1(glo(1)-ng_m1:,glo(2)-ng_m1:,glo(3)-ng_m1:)
    real(kind=dp_t), intent(in   ) :: m2(glo(1)-ng_m2:,glo(2)-ng_m2:,glo(3)-ng_m2:)
    real(kind=dp_t), intent(inout) :: inner_prod

    ! local
    integer :: i,j,k

    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             inner_prod = inner_prod + m1(i,j,k)*m2(i,j,k)
          end do
       end do
    end do

  end subroutine cc_inner_prod_3d
  
  subroutine stag_l2_norm(mla, m, norm)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:,:)     ! (1, 1:dm)
    real(kind=dp_t), intent(inout) :: norm
    
    ! local 
    real(kind=dp_t) :: inner_prod(1:mla%dim)
    integer :: dm

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stag_l2_norm")

    dm = mla%dim

    call stag_inner_prod(mla,m,1,m,1,inner_prod)
    norm = sqrt(sum(inner_prod(1:dm)))
    
    call destroy(bpt)

  end subroutine stag_l2_norm

  subroutine cc_l2_norm(mla, m, norm)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:)
    real(kind=dp_t), intent(inout) :: norm
    
    ! local 
    real(kind=dp_t) :: inner_prod

    type(bl_prof_timer), save :: bpt

    call build(bpt,"cc_l2_norm")

    call cc_inner_prod(mla,m,1,m,1,inner_prod)
    norm = sqrt(inner_prod)
    
    call destroy(bpt)

  end subroutine cc_l2_norm

  subroutine stag_l1_norm(mla,m,norm)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:,:)
    real(kind=dp_t), intent(inout) :: norm(:)

    ! local
    integer :: i,dm,nlevs,ng_m
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: norm_proc(mla%dim), norm_grid(mla%dim)

    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stag_l1_norm")

    norm_proc = 0.d0

    ng_m = m(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim

    if (nlevs .gt. 1) then
       call bl_error('staggered_norm not written for multilevel yet')
    end if

    do i=1,nfabs(m(1,1))
       mxp => dataptr(m(1,1), i)
       myp => dataptr(m(1,2), i)
       lo = lwb(get_box(m(1,1), i))
       hi = upb(get_box(m(1,1), i))
       norm_grid = 0.d0
       select case (dm)
       case (2)
          call stag_l1_norm_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                               lo, hi, norm_grid)
       case (3)
          mzp => dataptr(m(1,3), i)
          call stag_l1_norm_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                               lo, hi, norm_grid)
       end select

       norm_proc(1:dm) = norm_proc(1:dm) + norm_grid(1:dm)

    end do

    call parallel_reduce(norm(1:dm), norm_proc(1:dm), MPI_SUM)
    
    call destroy(bpt)

  end subroutine stag_l1_norm

  subroutine stag_l1_norm_2d(mx,my,ng_m,lo,hi,norm)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: norm(:)

    ! local
    integer :: i,j

    ! mx, interior cells
    do j=lo(2),hi(2)
       do i=lo(1)+1,hi(1)
          norm(1) = norm(1) + abs(mx(i,j))
       end do
    end do

    ! mx, boundary cells
    do j=lo(2),hi(2)
       norm(1) = norm(1) + 0.5d0*abs(mx(lo(1),j))
       norm(1) = norm(1) + 0.5d0*abs(mx(hi(1)+1,j))
    end do

    ! my, interior cells
    do j=lo(2)+1,hi(2)
       do i=lo(1),hi(1)
          norm(2) = norm(2) + abs(my(i,j))
       end do
    end do
    
    ! my, boundary cells
    do i=lo(1),hi(1)
       norm(2) = norm(2) + 0.5d0*abs(my(i,lo(2)))
       norm(2) = norm(2) + 0.5d0*abs(my(i,hi(2)+1))
    end do
    
     end subroutine stag_l1_norm_2d

  subroutine stag_l1_norm_3d(mx,my,mz,ng_m,lo,hi,norm)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) ::   mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: norm(:)

    ! local
    integer :: i,j,k

    ! mx, interior cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             norm(1) = norm(1) + abs(mx(i,j,k))
          end do
       end do
    end do

    ! mx, boundary cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          norm(1) = norm(1) + 0.5d0*abs(mx(lo(1),j,k))
          norm(1) = norm(1) + 0.5d0*abs(mx(hi(1)+1,j,k))
       end do
    end do

    ! my, interior cells
    do k=lo(3),hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             norm(2) = norm(2) + abs(my(i,j,k))
          end do
       end do
    end do

    ! my, boundary cells
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          norm(2) = norm(2) + 0.5d0*abs(my(i,lo(2),k))
          norm(2) = norm(2) + 0.5d0*abs(my(i,hi(2)+1,k))
       end do
    end do

    ! mz, interior cells
    do k=lo(3)+1,hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             norm(3) = norm(3) + abs(mz(i,j,k))
          end do
       end do
    end do

    ! mz, boundary cells
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          norm(3) = norm(3) + 0.5d0*abs(mz(i,j,lo(3)))
          norm(3) = norm(3) + 0.5d0*abs(mz(i,j,hi(3)+1))
       end do
    end do

  end subroutine stag_l1_norm_3d

  subroutine sum_umac_press(mla,pressure,umac,av_sump,av_sumu)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: pressure(:)
    type(multifab) , intent(in   ) ::    umac(:,:)
    real(kind=dp_t), intent(out), optional :: av_sump, av_sumu(mla%dim) ! Level 1 only

    ! local
    integer :: i,dm,nlevs,ng_m
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: sump_lev
    real(kind=dp_t) :: sumu_lev(mla%dim), sumu_proc(mla%dim), sumu_grid(mla%dim)

    real(kind=dp_t), pointer ::  cp(:,:,:,:)
    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"sum_umac_press")

    sumu_lev = 0.d0
    sumu_proc = 0.d0

    ng_m = umac(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim
    
    if (nlevs .gt. 1) then
       call bl_error('sum_umac_press not written for multilevel yet')
    end if

    ! add up pressure
    sump_lev = multifab_sum_c(pressure(1),1,1,all=.false.)

    !$omp parallel private(mfi,i,xnodalbox,xlo,xhi,ynodalbox,ylo,yhi) &
    !$omp private(znodalbox,zlo,zhi,cp,mxp,myp,lo,hi,sumu_grid) &
    !$omp reduction(+:sumu_proc)

    ! add up umac
    call mfiter_build(mfi, pressure(1), tiling=.true.)

    do while (more_tile(mfi))
       i = get_fab_index(mfi)
       
       xnodalbox = get_nodaltilebox(mfi,1)
       xlo = lwb(xnodalbox)
       xhi = upb(xnodalbox)
       ynodalbox = get_nodaltilebox(mfi,2)
       ylo = lwb(ynodalbox)
       yhi = upb(ynodalbox)
       znodalbox = get_nodaltilebox(mfi,3)
       zlo = lwb(znodalbox)
       zhi = upb(znodalbox)

!    do i=1,nfabs(pressure(1))
       cp  => dataptr(pressure(1), i)
       mxp => dataptr(umac(1,1), i)
       myp => dataptr(umac(1,2), i)
       lo = lwb(get_box(pressure(1), i))
       hi = upb(get_box(pressure(1), i))
       sumu_grid = 0.d0
       select case (dm)
       case (2)
          call sum_umac_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                           lo, hi, sumu_grid,xlo,xhi,ylo,yhi)
       case (3)
          mzp => dataptr(umac(1,3), i)
          call sum_umac_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                           lo, hi, sumu_grid,xlo,xhi,ylo,yhi,zlo,zhi)

       end select

       sumu_proc(1:dm) = sumu_proc(1:dm) + sumu_grid(1:dm)

    end do
    !$omp end parallel

    call parallel_reduce(sumu_lev(1:dm), sumu_proc(1:dm), MPI_SUM)

    if (parallel_IOProcessor().and.(.not.present(av_sump))) then
       write(*,"(A,100G17.9)") "<p>=", sump_lev/total_volume
    end if
    if (parallel_IOProcessor().and.(.not.present(av_sumu))) then       
       write(*,"(A,100G17.9)") "<u>=", sumu_lev(1:dm)/total_volume
    end if
        
    if(present(av_sump)) av_sump=sump_lev/total_volume
    if(present(av_sumu)) av_sumu(1:dm)=sumu_lev(1:dm)/total_volume
    
    call destroy(bpt)

  end subroutine sum_umac_press

  subroutine sum_umac_2d(umac,vmac,ng_m,glo,ghi,sumu,xlo,xhi,ylo,yhi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m
    integer        , intent(in   ) :: xlo(:), xhi(:), ylo(:), yhi(:)
    real(kind=dp_t), intent(in   ) ::   umac(glo(1)-ng_m:,glo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   vmac(glo(1)-ng_m:,glo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: sumu(:)

    ! local
    integer :: i,j

    ! umac, interior cells
    do j=xlo(2),xhi(2)
       do i=xlo(1),xhi(1)
          sumu(1) = sumu(1) + umac(i,j)
       end do
    end do

    ! umac, boundary cells
    if(xlo(1) .eq. glo(1) .and. xhi(1) .eq. ghi(1)+1) then
    do j=xlo(2),xhi(2)
       sumu(1) = sumu(1) + 0.5d0*umac(xlo(1)  ,j)
       sumu(1) = sumu(1) + 0.5d0*umac(xhi(1),j)
    end do
    end if

    ! vmac, interior cells
    do j=ylo(2),yhi(2)
       do i=ylo(1),yhi(1)
          sumu(2) = sumu(2) + vmac(i,j)
       end do
    end do

    ! vmac, boundary cells
    if(ylo(2) .eq. glo(2) .and. yhi(2) .eq. ghi(2)+1) then
    do i=ylo(1),yhi(1)
       sumu(2) = sumu(2) + 0.5d0*vmac(i,ylo(2)  )
       sumu(2) = sumu(2) + 0.5d0*vmac(i,yhi(2))
    end do
    end if

  end subroutine sum_umac_2d

  subroutine sum_umac_3d(umac,vmac,wmac,ng_m,glo,ghi,sumu,xlo,xhi,ylo,yhi,zlo,zhi)

    integer        , intent(in   ) :: glo(:), ghi(:), ng_m
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t), intent(in   ) ::   umac(glo(1)-ng_m:,glo(2)-ng_m:,glo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   vmac(glo(1)-ng_m:,glo(2)-ng_m:,glo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) ::   wmac(glo(1)-ng_m:,glo(2)-ng_m:,glo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: sumu(:)

    ! local
    integer :: i,j,k

    ! umac, interior cells
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)
             sumu(1) = sumu(1) + umac(i,j,k)
          end do
       end do
    end do

    ! umac, boundary cells
    if(xlo(1) .eq. glo(1) .and. xhi(1) .eq. ghi(1)+1) then
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          sumu(1) = sumu(1) + 0.5d0*umac(xlo(1)  ,j,k)
          sumu(1) = sumu(1) + 0.5d0*umac(xhi(1),j,k)
       end do
    end do
    end if

    ! vmac, interior cells
    do k=ylo(3),yhi(3)
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)
             sumu(2) = sumu(2) + vmac(i,j,k)
          end do
       end do
    end do

    ! vmac, boundary cells
    if(ylo(2) .eq. glo(2) .and. yhi(2) .eq. ghi(2)+1) then
    do k=ylo(3),yhi(3)
       do i=ylo(1),yhi(1)
          sumu(2) = sumu(2) + 0.5d0*vmac(i,ylo(2)  ,k)
          sumu(2) = sumu(2) + 0.5d0*vmac(i,yhi(2),k)
       end do
    end do
    end if

    ! wmac, interior cells
    do k=zlo(3),zhi(3)
       do j=zlo(2),zhi(2)
          do i=zlo(1),zhi(1)
             sumu(3) = sumu(3) + wmac(i,j,k)
          end do
       end do
    end do

    ! wmac, boundary cells
    if(zlo(3) .eq. glo(3) .and. zhi(3) .eq. ghi(3)+1) then
    do j=zlo(2),zhi(2)
       do i=zlo(1),zhi(1)
          sumu(3) = sumu(3) + 0.5d0*wmac(i,j,zlo(3)  )
          sumu(3) = sumu(3) + 0.5d0*wmac(i,j,zhi(3))
       end do
    end do
    end if

  end subroutine sum_umac_3d

end module norm_inner_product_module
