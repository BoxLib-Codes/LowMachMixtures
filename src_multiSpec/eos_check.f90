module eos_check_module

  use multifab_module
  use ml_layout_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: rhobar

  implicit none

  private

  public :: eos_check

contains

  subroutine eos_check(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)

    ! local
    integer i,n,dm,nlevs,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: rp(:,:,:,:)

    real(kind=dp_t) :: eos_error, eos_error_grid, eos_error_proc

    type(bl_prof_timer), save :: bpt

    call build(bpt, "eos_check")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_r = rho(1)%ng

    eos_error_proc = -1.d20
    do n=1,nlevs
       do i=1,nfabs(rho(n))
         rp => dataptr(rho(n), i)
         lo =  lwb(get_box(rho(n), i))
         hi =  upb(get_box(rho(n), i))
         eos_error_grid = -1.d20
         select case (dm)
         case (2)
            call eos_check_2d(rp(:,:,1,:),ng_r,eos_error_grid,lo,hi)
         case (3)
            call eos_check_3d(rp(:,:,:,:),ng_r,eos_error_grid,lo,hi)
         end select
         eos_error_proc = max(eos_error_grid, eos_error_proc)
      end do
   end do

   ! This sets eos_error to be the max of eos_error_proc over all processors.
   call parallel_reduce(eos_error, eos_error_proc, MPI_MAX)

   if (parallel_IOProcessor()) then
      print*,"EOS ERROR in L1 norm: ",eos_error
      print*,""
   end if

   call destroy(bpt)

  end subroutine eos_check

  subroutine eos_check_2d(rho,ng_r,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,n

    real(kind=dp_t) :: sum,error

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,n)/rhobar(n)
          end do
          error = abs(1.d0 - sum)
          eos_error = max(eos_error,error)

       end do
    end do

  end subroutine eos_check_2d

  subroutine eos_check_3d(rho,ng_r,eos_error,lo,hi)

    integer        , intent(in   ) :: lo(:), hi(:), ng_r
    real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t), intent(inout) :: eos_error

    ! local
    integer :: i,j,k,n

    real(kind=dp_t) :: sum,error

    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             sum = sum + rho(i,j,k,n)/rhobar(n)
          end do
          error = abs(1.d0 - sum)
          eos_error = max(eos_error,error)

             
          end do
       end do
    end do

  end subroutine eos_check_3d

end module eos_check_module
