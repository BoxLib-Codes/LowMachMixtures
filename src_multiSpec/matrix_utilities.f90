module matrix_utilities

  use bl_types
  use bl_prof_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: molmass

  implicit none

  private

  public :: Dbar2chi_iterative, choldc
  
contains

    ! nspecies is number of species
    ! num_iterations is the number of terms in the sum to use: 3-5 are reasonable values
    ! D_bar is matrix of Maxwell-Stefan binary diffusion coefficient
    ! chi is the multispecies diffusion matrix
    ! Xk is mole fractions --- MUST NOT BE ZERO
    subroutine Dbar2chi_iterative(num_iterations,D_bar,Xk,chi)
      integer, intent(in) :: num_iterations
      real(kind=dp_t), intent(in) :: Xk(1:nspecies), D_bar(1:nspecies,1:nspecies)
      real(kind=dp_t), intent(out) :: chi(1:nspecies,1:nspecies)
      
      ! Local variables
      real(kind=dp_t) :: term1, term2, MWmix
      real(kind=dp_t) :: Di(1:nspecies)
      real(kind=dp_t) :: Deltamat(1:nspecies,1:nspecies), Zmat(1:nspecies,1:nspecies)
      real(kind=dp_t), dimension(1:nspecies,1:nspecies) :: Pmat, Jmat
      real(kind=dp_t), dimension(1:nspecies) :: Minv, Mmat
      real(kind=dp_t), dimension(1:nspecies,1:nspecies) :: PJ, matrix1, matrix2
      real(kind=dp_t) :: scr
      real(kind=dp_t) :: Ykp(1:nspecies), Xkp(1:nspecies)

      integer :: i, j, k, ii, jj

      type(bl_prof_timer), save :: bpt

      call build(bpt,"Dbar2chi_iterative")

      ! mole fractions correction
      ! Turned this off since it should be done in the caller
      do ii = 1, nspecies
       Xkp(ii) = Xk(ii) ! + fraction_tolerance*(sum(Xk(:))/dble(nspecies)-Xk(ii))
      end do

      ! molecular weight of mixture - EGLIB
      Mwmix = 0.0d0
      do ii = 1, nspecies
       MWmix = MWmix + Xkp(ii)*molmass(ii)
      end do

      ! mass fractions correction - EGLIB
      do ii = 1, nspecies
       Ykp(ii) = molmass(ii)/MWmix*Xkp(ii)
      end do

      ! Find Di matrix 
      do i = 1, nspecies
       term2 = 0.0d0
       do j = 1, nspecies
        if(j.ne.i) then
          term2 = term2 + Xkp(j)/D_bar(i,j)
        end if
       end do   
       Di(i) = (1.d0-Ykp(i))/term2 
      end do   

      ! Compute Mmat and Minv
      do i = 1, nspecies
       Mmat(i) = Xkp(i)/Di(i)
       Minv(i) = Di(i)/Xkp(i)
      end do

      ! Compute P matrix
      Pmat = 0.0d0
      do i = 1, nspecies
       do j = 1, nspecies
         Pmat(i,j) = - Ykp(j) 
         if(i.eq.j) then
          Pmat(i,j) =  Pmat(i,j) + 1.0d0  
         end if
       end do
      end do

      ! Compute Deltamat
      Deltamat = 0.0d0 
      do i = 1, nspecies
       do j = 1, nspecies
         if(i.eq.j) then
          term1 = 0.0d0
          do k = 1, nspecies
           if(k.ne.i) then
            term1 = term1 + Xkp(i)*Xkp(k)/D_bar(i,k)
           end if
          end do  
          Deltamat(i,i) = term1
         else
          Deltamat(i,j) = -Xkp(i)*Xkp(j)/D_bar(i,j) 
         end if  
          Zmat(i,j) = -Deltamat(i,j)
       end do
      end do  

      ! Compute Zmat
      do i = 1, nspecies
        Zmat(i,i) = Zmat(i,i) + Mmat(i)
      end do  

      ! Compute Jmat
      do i = 1, nspecies
       do j = 1, nspecies
         Jmat(i,j) = Minv(i)*Zmat(i,j)
        end do
       end do

      ! Compute PJ
      PJ = 0.0d0
      do i = 1, nspecies
       do j = 1, nspecies
        do k = 1, nspecies
         PJ(i,j) = PJ(i,j) + Pmat(i,k)*Jmat(k,j)
        end do
       end do
      end do

      ! Compute P M^-1 Pt; store it in matrix2
      do i = 1, nspecies
       do j = 1, nspecies
        scr = 0.d0
        do k = 1, nspecies
         scr = scr + Pmat(i,k)*Minv(k)*Pmat(j,k) 
            ! notice the change in indices for Pmat to represent Pmat^t
        end do
         matrix2(i,j) = scr
         chi(i,j) = scr
       end do
      end do

      do jj = 1,num_iterations
       do i = 1, nspecies
        do j = 1, nspecies
         scr = 0.d0
         do k = 1, nspecies
            scr = scr + PJ(i,k)*chi(k,j)
         end do
          matrix1(i,j) = scr+matrix2(i,j)
        end do
       end do 
       chi=matrix1
      end do

      call destroy(bpt)

  end subroutine

   ! a is input matrix.  
   ! upon return the lower triangle and diagonal are overwritten by the cholesky factor
   subroutine choldc(a,np)
       integer :: np
       real(kind=dp_t), intent(inout) :: a(np,np)

       real(kind=dp_t) :: p(np), dij(np,np)
       real(kind=dp_t) :: dd(np,np)
       real(kind=dp_t) :: yy(np), mwmix 

       integer :: i, j, k, ii, jj
       real(kind=dp_t) :: sum1
       real(kind=dp_t) :: small_number = 0.0d0 ! Some tolerance

       integer :: idiag,ising

       type(bl_prof_timer), save :: bpt

       call build(bpt,"choldc")

       do i = 1, np

           ising = 0

        do j = i, np

           sum1 = a(i,j)

           do k = i-1, 1, -1

              sum1 = sum1 - a(i,k)*a(j,k)

           end do

           if(i.eq.j) then

             if(sum1.le.small_number) then

             p(i) = 0.d0

             ising = 1

             else

             p(i) = sqrt(sum1)

             end if

           else

             if(ising.eq.0)then

                a(j,i) = sum1/p(i)

             else

                a(j,i) = 0.d0

             end if

           end if

        end do

       end do


       do i = 1, np

          do j = i+1, np

           a(i,j) = 0.0d0 ! Zero upper triangle

          end do
          
          a(i,i) = p(i)

       end do

       call destroy(bpt)

    end subroutine

end module
