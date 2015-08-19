module sum_momenta_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: total_volume

  implicit none

  private

  public :: sum_momenta

contains

  subroutine sum_momenta(mla,m,av_m)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: m(:,:)
    real(kind=dp_t), intent(inout), optional :: av_m(:)

    ! local
    integer :: i,dm,nlevs,ng_m
    integer :: lo(mla%dim), hi(mla%dim)

    real(kind=dp_t) :: mom_tot(mla%dim), mom_lev(mla%dim), mom_proc(mla%dim), mom_grid(mla%dim)

    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"sum_momenta")

    mom_tot = 0.d0
    mom_lev = 0.d0
    mom_proc = 0.d0

    ng_m = m(1,1)%ng
    
    nlevs = mla%nlevel
    dm = mla%dim
    
    if (nlevs .gt. 1) then
       call bl_error('sum_momenta not written for multilevel yet')
    end if

    do i=1,nfabs(m(1,1))
       mxp => dataptr(m(1,1), i)
       myp => dataptr(m(1,2), i)
       lo = lwb(get_box(m(1,1), i))
       hi = upb(get_box(m(1,1), i))
       mom_grid  = 0.d0
       select case (dm)
       case (2)
          call sum_momenta_2d(mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                              lo, hi, mom_grid)
       case (3)
          mzp => dataptr(m(1,3), i)
          call sum_momenta_3d(mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                              lo, hi, mom_grid)

       end select

       mom_proc(1:dm) = mom_proc(1:dm) + mom_grid(1:dm)

    end do

    call parallel_reduce(mom_lev(1:dm)    , mom_proc(1:dm)    , MPI_SUM)

    if (parallel_IOProcessor()) then
       write(*,"(A,100G17.9)") "CONSERVE: <mom_k>=", mom_lev(1:dm)/total_volume
    end if
        
    if(present(av_m)) av_m=mom_lev(1:dm)/total_volume
    
    call destroy(bpt)

  end subroutine sum_momenta

  subroutine sum_momenta_2d(mx,my,ng_m,lo,hi,mom)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) :: mx(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(in   ) :: my(lo(1)-ng_m:,lo(2)-ng_m:)
    real(kind=dp_t), intent(inout) :: mom(:)

    ! local
    integer :: i,j

    ! mx, interior cells
    do j=lo(2),hi(2)
       do i=lo(1)+1,hi(1)
          mom(1) = mom(1) + mx(i,j)
       end do
    end do

    ! mx, boundary cells
    do j=lo(2),hi(2)
       mom(1) = mom(1) + 0.5d0*mx(lo(1)  ,j)
       mom(1) = mom(1) + 0.5d0*mx(hi(1)+1,j)
    end do

    ! my, interior cells
    do j=lo(2)+1,hi(2)
       do i=lo(1),hi(1)
          mom(2) = mom(2) + my(i,j)
       end do
    end do

    ! my, boundary cells
    do i=lo(1),hi(1)
       mom(2) = mom(2) + 0.5d0*my(i,lo(2)  )
       mom(2) = mom(2) + 0.5d0*my(i,hi(2)+1)
    end do

  end subroutine sum_momenta_2d

  subroutine sum_momenta_3d(mx,my,mz,ng_m,lo,hi,mom)

    integer        , intent(in   ) :: lo(:), hi(:), ng_m
    real(kind=dp_t), intent(in   ) :: mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) :: my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(in   ) :: mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
    real(kind=dp_t), intent(inout) :: mom(:)

    ! local
    integer :: i,j,k

    ! mx, interior cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1)+1,hi(1)
             mom(1) = mom(1) + mx(i,j,k)
          end do
       end do
    end do

    ! mx, boundary cells
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          mom(1) = mom(1) + 0.5d0*mx(lo(1)  ,j,k)
          mom(1) = mom(1) + 0.5d0*mx(hi(1)+1,j,k)
       end do
    end do

    ! my, interior cells
    do k=lo(3),hi(3)
       do j=lo(2)+1,hi(2)
          do i=lo(1),hi(1)
             mom(2) = mom(2) + my(i,j,k)
          end do
       end do
    end do

    ! my, boundary cells
    do k=lo(3),hi(3)
       do i=lo(1),hi(1)
          mom(2) = mom(2) + 0.5d0*my(i,lo(2)  ,k)
          mom(2) = mom(2) + 0.5d0*my(i,hi(2)+1,k)
       end do
    end do

    ! mz, interior cells
    do k=lo(3)+1,hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             mom(3) = mom(3) + mz(i,j,k)
          end do
       end do
    end do

    ! mz, boundary cells
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          mom(3) = mom(3) + 0.5d0*mz(i,j,lo(3)  )
          mom(3) = mom(3) + 0.5d0*mz(i,j,hi(3)+1)
       end do
    end do

  end subroutine sum_momenta_3d
  
end module sum_momenta_module
