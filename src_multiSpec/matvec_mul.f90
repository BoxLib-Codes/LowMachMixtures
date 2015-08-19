module matvec_mul_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: matvec_mul

  contains

  subroutine matvec_mul(mla, x, A, nc) ! Computes x=A*x, where A is an nxn matrix and x is an nx1 vector  
   
    type(ml_layout), intent(in  )  :: mla
    type(multifab),  intent(inout) :: x
    type(multifab),  intent(in)    :: A
    integer,         intent(in)    :: nc
 
    ! local variables
    real(kind=dp_t), pointer       :: xp(:,:,:,:) 
    real(kind=dp_t), pointer       :: ap(:,:,:,:)
    integer                        :: lo(mla%dim), hi(mla%dim), i, dm
   
    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"matvec_mul")

    dm = mla%dim 

    !$omp parallel private(mfi,i,tilebox,tlo,thi,xp,ap,lo,hi)

    call mfiter_build(mfi, A, tiling=.true.)

    do while(more_tile(mfi))
       i = get_fab_index(mfi)

       tilebox = get_tilebox(mfi)
       tlo = lwb(tilebox)
       thi = upb(tilebox)

!    do i=1,nfabs(x)
       xp  => dataptr(x, i)   
       ap  => dataptr(A, i)  
       lo(1) = lbound(ap,1)  ! this adjusts lo & hi ghost cells on nodes & faces,
       lo(2) = lbound(ap,2)  ! so in the subroutine loop from lo to hi without ng.
       hi(1) = ubound(ap,1)
       hi(2) = ubound(ap,2)

       select case (dm)
         case (2)
             call matvec_mul_2d(xp(:,:,1,:), ap(:,:,1,:), lo, hi, nc, tlo, thi)
         case (3)
             lo(3) = lbound(ap,3) 
             hi(3) = ubound(ap,3)
             call matvec_mul_3d(xp(:,:,:,:), ap(:,:,:,:), lo, hi, nc, tlo, thi)
       end select
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine matvec_mul

  subroutine matvec_mul_2d(xp, ap, lo, hi, nc, tlo, thi)

    integer                        :: lo(:), hi(:), tlo(:), thi(:)
    real(kind=dp_t), intent(inout) :: xp(lo(1):,lo(2):,:) ! last dimension for nc
    real(kind=dp_t), intent(in)    :: ap(lo(1):,lo(2):,:)
    integer                        :: i,j,nc

    !print*, lo(1:2), hi(1:2)
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
          call matvec_mul_comp(xp(i,j,:), ap(i,j,:))
          !if(i.eq.7 .and. j.eq.14) print*, 'x-flux / y-flux'
          !if(i.eq.7 .and. j.eq.14) print*, "flux1=",xp(i,j,1), "flux2=",xp(i,j,2),"flux1+flux2=",xp(i,j,1)+xp(i,j,2) 
          !print*, i,j,"x-flux=",xp(i,j,1)!, "x-flux=",xp(i,j,2)
       end do
    end do
          !print*, '' 
          !print*, "Now y-flux for 2-species"

    contains 
    
    ! Use contained subroutine to do rank conversion and mat-vec mult
    subroutine matvec_mul_comp(xp_ij, ap_ij)

        real(kind=dp_t), dimension(nc),    intent(inout) :: xp_ij
        real(kind=dp_t), dimension(nc,nc), intent(in)    :: ap_ij  
        
        xp_ij = matmul(ap_ij, xp_ij)
 
    end subroutine matvec_mul_comp 

  end subroutine matvec_mul_2d
 
  subroutine matvec_mul_3d(xp, ap, lo, hi, nc, tlo, thi)

    integer                        :: lo(:), hi(:), tlo(:), thi(:)
    real(kind=dp_t), intent(inout) :: xp(lo(1):,lo(2):,lo(3):,:) ! last dimension for nc
    real(kind=dp_t), intent(in)    :: ap(lo(1):,lo(2):,lo(3):,:)
    integer                        :: i,j,k,nc

    !print*, lo(1:3), hi(1:3)
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
             call matvec_mul_comp(xp(i,j,k,:), ap(i,j,k,:))
             !if(i.eq.7 .and. j.eq.14) print*, 'x-flux / y-flux'
             !if(i.eq.7 .and. j.eq.14) print*, "flux1=",xp(i,j,1), "flux2=",xp(i,j,2),"flux1+flux2=",xp(i,j,1)+xp(i,j,2) 
             !print*, i,j,"x-flux=",xp(i,j,1)!, "x-flux=",xp(i,j,2)
          end do
       end do
    end do
             !print*, '' 
             !print*, "Now y-flux for 2-species"
    contains 
    
    ! Use contained subroutine to do rank conversion and mat-vec mult
    subroutine matvec_mul_comp(xp_ij, ap_ij)

        real(kind=dp_t), dimension(nc),    intent(inout) :: xp_ij
        real(kind=dp_t), dimension(nc,nc), intent(in)    :: ap_ij  
        
        xp_ij = matmul(ap_ij, xp_ij)
 
    end subroutine matvec_mul_comp 

  end subroutine matvec_mul_3d

end module matvec_mul_module
