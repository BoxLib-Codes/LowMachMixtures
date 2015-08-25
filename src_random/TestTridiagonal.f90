program TestTridiagonal

   use, intrinsic :: ISO_C_BINDING ! For interfacing with C
   use Precision
   use TridiagonalSolvers

   integer, parameter :: wp=c_double ! Double precision

   integer, parameter :: n=10
   
   integer :: i, j, k
   real(wp) :: matrix(n,n), u(0:n+1)
   real(wp), dimension(n) :: a, b, c, r, residual
   
   a=1.0_wp
   b=-2.0_wp
   c=1.0_wp
   
   do i=1,n     
      r(i)=i
   end do
      
   call SolveTridiagonal(a,b,c,r,u(1:n),n,periodic=0)

   u(0)=0
   u(n+1)=0
   do i=1,n
      residual(i) = (u(i-1)-2*u(i)+u(i+1)) - r(i)
   end do
   
   write(*,*) "STENCIL Residual = ", maxval(abs(residual))

   do i=1,n     
      a(i)=2*n+i
      b(i)=i
      c(i)=2*n-i      
   end do

   call SolveTridiagonal(a,b,c,r,u(1:n),n,periodic=0)
   
   do i=1,n
      if(i>1) matrix(i,i-1)=a(i)
      matrix(i,i)=b(i)
      if(i>1) matrix(i-1,i)=c(i-1)      
   end do
   
   residual = matmul(matrix, u(1:n)) - r

   write(*,*) "MATVEC Residual = ", maxval(abs(residual))
   
end program
