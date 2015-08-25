! *************************************************************************
! ** Tridiagonal (1D) solver from John Bell 
! ** Do a tridiagonal solve Au=r, where A=[a,b,c] with periodic or (zero) Dirichlet BCs
! The C interface for the routines here is in HydroGrid.h
!
! *************************************************************************
module TridiagonalSolvers
   use, intrinsic :: ISO_C_BINDING ! For interfacing with C
   use Precision   
   implicit none
   private
   public :: SolveTridiagonal
   
   integer, parameter :: wp=c_double ! Double precision

contains

! The C prototype for this routine is
! void SolveTridiagonal(double *a, double *b, double *c, double *r, double *u, int n, int periodic)
subroutine SolveTridiagonal(a,b,c,r,u,n,periodic) bind(C, name="SolveTridiagonal")
      implicit none      
      integer(c_int), value :: n
      real(wp), intent(in), dimension(n) :: a, b, c, r
      real(wp), intent(out), dimension(n) :: u
      integer(c_int), value :: periodic

      real(wp), dimension(n) :: ll, ur, air, aib, gam
      real(wp) bet, alpha, delta

      integer j
      integer nmax

      if (b(1) .eq. 0) then
         print *,'CANT HAVE B(1) = ZERO'
         stop
      end if

      do j=1,n-1
         ur(j) = 0.0_wp
         ll(j) = 0.0_wp
      enddo

      ur(n-1) = c(n-1)
      
      if(periodic>0) then
         ur(1) = a(1)
         ll(1) = c(n)
         ll(n-1) = a(n)
      else
         ll(n-1) = a(n)
      end if
      
      bet = b(1)
      air(1) = r(1)/bet
      aib(1) = ur(1)/bet

      do j = 2,n-1
        gam(j) = c(j-1)/bet
        bet = b(j) - a(j)*gam(j)
        if (bet .eq. 0) then
          print *,'TRIDIAG FAILED '
          stop
        endif
        air(j) = (r(j)-a(j)*air(j-1))/bet
        aib(j) = (ur(j)-a(j)*aib(j-1))/bet
      enddo

      do j = n-2,1,-1
        air(j) = air(j) - gam(j+1)*air(j+1)
        aib(j) = aib(j) - gam(j+1)*aib(j+1)
      enddo

      alpha = 0.0_wp
      delta = 0.0_wp

      do j=1,n-1   
          alpha = alpha + ll(j)*air(j)
          delta = delta + ll(j)*aib(j)
      enddo

      u(n) = (r(n) - alpha) / (b(n)-delta)
      
      do j=1,n-1
         u(j) = air(j) - u(n)*aib(j)
      enddo

end subroutine
        
end module
