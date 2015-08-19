module matmat_mul_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: matmat_mul

  contains

  subroutine matmat_mul(mla, x, A, nc) ! Computes x=A*x, where A and x are ncXnc matrix 
   
    type(ml_layout), intent(in   ) :: mla
    type(multifab),  intent(inout) :: x
    type(multifab),  intent(in   ) :: A
    integer,         intent(in   ) :: nc
 
    ! local variables
    real(kind=dp_t), pointer       :: xp(:,:,:,:) 
    real(kind=dp_t), pointer       :: ap(:,:,:,:)
    integer                        :: lo(mla%dim), hi(mla%dim), i, dm
   
    type(bl_prof_timer), save :: bpt

    call build(bpt,"matmat_mul")

    dm = mla%dim 
 
    do i=1,nfabs(x)
       xp  => dataptr(x, i)   
       ap  => dataptr(A, i)  
       lo(1) = lbound(ap,1)  ! this adjusts lo & hi ghost cells on nodes & faces,
       lo(2) = lbound(ap,2)  ! so in the subroutine loop from lo to hi without ng.
       hi(1) = ubound(ap,1)
       hi(2) = ubound(ap,2)
       !print*, 'matmat', lo(1:2),hi(1:2)

       select case (dm)
         case (2)
             call matmat_mul_2d(xp(:,:,1,:), ap(:,:,1,:), lo, hi, nc)
         case (3)
             lo(3) = lbound(ap,3) 
             hi(3) = ubound(ap,3)
             call matmat_mul_3d(xp(:,:,:,:), ap(:,:,:,:), lo, hi, nc)
       end select
    end do

    call destroy(bpt)

  end subroutine matmat_mul

  subroutine matmat_mul_2d(xp, ap, lo, hi, nc)

    integer                        :: lo(:), hi(:) 
    real(kind=dp_t), intent(inout) :: xp(lo(1):,lo(2):,:) ! last dimension for nc^2
    real(kind=dp_t), intent(in)    :: ap(lo(1):,lo(2):,:)
    integer                        :: i,j,nc

    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
          call matmat_mul_comp(xp(i,j,:), ap(i,j,:))
       end do
    end do

    contains 
    
    ! Use contained subroutine to do rank conversion and matrix-matrix multiplication
    subroutine matmat_mul_comp(xp_ij, ap_ij)

        real(kind=dp_t), dimension(nc,nc), intent(inout) :: xp_ij
        real(kind=dp_t), dimension(nc,nc), intent(in)    :: ap_ij  
        !print*, nc 
        xp_ij = matmul(ap_ij, xp_ij)
  
    end subroutine matmat_mul_comp 

  end subroutine matmat_mul_2d
 
  subroutine matmat_mul_3d(xp, ap, lo, hi, nc)

    integer                        :: lo(:), hi(:) 
    real(kind=dp_t), intent(inout) :: xp(lo(1):,lo(2):,lo(3):,:) ! last dimension for nc^2
    real(kind=dp_t), intent(in)    :: ap(lo(1):,lo(2):,lo(3):,:)
    integer                        :: i,j,k,nc

    !print*, lo(1:3), hi(1:3)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)
             call matmat_mul_comp(xp(i,j,k,:), ap(i,j,k,:))
          end do
       end do
    end do
    
    contains 
    
    ! Use contained subroutine to do rank conversion and matrix-matrix multiplication
    subroutine matmat_mul_comp(xp_ij, ap_ij)

        real(kind=dp_t), dimension(nc,nc), intent(inout) :: xp_ij
        real(kind=dp_t), dimension(nc,nc), intent(in)    :: ap_ij  
        
        type(bl_prof_timer), save :: bpt

        call build(bpt,"matmat_mul_comp")

        xp_ij = matmul(ap_ij, xp_ij)
 
        call destroy(bpt)

    end subroutine matmat_mul_comp 

  end subroutine matmat_mul_3d

end module matmat_mul_module
