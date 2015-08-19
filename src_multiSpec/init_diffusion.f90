module init_diffusion_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use multifab_coefbc_module
  use ml_layout_module
  use convert_stag_module
  use convert_rhoc_to_c_module
  use probin_common_module, only: prob_lo, prob_hi, prob_type, molmass
  use probin_multispecies_module, only: alpha1, beta, delta, sigma, Dbar, &
                                        c_init, nspecies
 
  implicit none

  private

  public :: init_rho  ! used in diffusion code

  ! IMPORTANT: In the diffusion only code (init_rho), c_init specifies initial values for DENSITY
  ! In the low-Mach code (init_rho_and_umac), c_init specifies initial MASS FRACTIONS (should sum to unity!)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_type codes for LowMach:
  ! 0=thermodynamic equilibrium, v=0, rho/c=c_init(1,1:nspecies)
  ! 1=bubble test, v=0, rho/c=c_init(1,1:nspecies) inside and c_init(2,1:nspecies) outside the bubble
  ! 2=gradient along y, v=0, rho/c=c_init(1,1:nspecies) on bottom (y=0) and c_init(2,1:nspecies) on top (y=Ly)
  ! 3=one fluid on top of another: v=0, rho/c=c_init(1,1:nspecies) on bottom (y<Ly/2) and c_init(2,1:nspecies) on top (y=L_y)
  ! 4=reserved for future use generic case (any nspecies)
  ! 5-onward=manufactured solutions for testing, limited to specific setups and nspecies

  
contains

  subroutine init_rho(rho,dx,time,the_bc_level)

    type(multifab) , intent(inout) :: rho(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: time
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer                        :: dm, ng_r, i, n, nlevs
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! pointer for rho (last dim:nspecies)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_rho")

    dm = rho(1)%dim
    ng_r = rho(1)%ng
    nlevs = size(rho,1)

    ! assign values of parameters for the gaussian rho, rhototal
    alpha1 = 0.5d0 
    beta   = 0.1d0 
    delta  = 0.5d0 
    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call init_rho_2d(dp(:,:,1,:),ng_r,lo,hi,dx(n,:),time)
          case (3)
             call init_rho_3d(dp(:,:,:,:),ng_r,lo,hi,dx(n,:),time)
          end select
       end do

       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(rho(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(rho(n),1,c_bc_comp,nspecies,the_bc_level(n), &
                            dx_in=dx(n,:))

    end do

    call destroy(bpt)

  end subroutine init_rho

  subroutine init_rho_2d(rho,ng_r,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng_r
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)  ! last dimension for species
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j,n
    real(kind=dp_t)  :: x,y,w1,w2,rsq,rhot,L(2),sum,r,y1,rho_loc,rand
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    select case (prob_type)

    case(0) 
    !============================================================
    ! Thermodynamic equilibrium
    !============================================================
 
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)
 
          rho(i,j,1:nspecies) = c_init(1,1:nspecies)

       end do
    end do  
    
    case(1) 
    !=============================================================
    ! Initializing rho's in concentric circle with radius^2 = 0.1
    !=============================================================
 
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
          rsq = x**2 + y**2
          if (rsq .lt. L(1)*L(2)*0.1d0) then
             rho(i,j,1:nspecies) = c_init(1,1:nspecies)
          else
             rho(i,j,1:nspecies) = c_init(2,1:nspecies)
          end if
    
       end do
    end do
  
    case(2) 
    !=========================================================
    ! Initializing rho's with constant gradient 
    !=========================================================
 
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
   
            ! linear gradient in rho
            rho(i,j,1:nspecies) = c_init(1,1:nspecies) + & 
               (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)
   
         end do
      end do

    case(3) 
    !===========================================================
    ! Initializing rho's in Gaussian so as rhotot=constant=1.0
    ! Here rho_exact = e^(-r^2/4Dt)/(4piDt)
    !===========================================================
 
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half*(prob_lo(2)+prob_hi(2))
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half*(prob_lo(1)+prob_hi(1))
        
            rsq = x**2 + y**2
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time))
            rho(i,j,2) = 1.0d0-1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time))
       
         end do
      end do

    case(4)
    !==================================================================================
    ! Initializing rho1,rho2=Gaussian and rhototal=1+alpha*exp(-r^2/4D)/(4piD) (no-time 
    ! dependence). Manufactured solution rho1_exact = exp(-r^2/4Dt-beta*t)/(4piDt)
    !==================================================================================
 
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2
            rhot = 1.0d0 + alpha1/(4.0d0*M_PI*Dbar(1))*dexp(-rsq/(4.0d0*Dbar(1)))
            rho(i,j,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)*dexp(-rsq/(4.0d0*Dbar(1)*time)-&
                         beta*time)*rhot
            rho(i,j,2) = rhot - rho(i,j,1)

         end do
      end do

    case(5)
    !==================================================================================
    ! Initializing m2=m3, D12=D13 where Dbar(1)=D12, Dbar(2)=D13, 
    ! Dbar(3)=D23, Grad(w2)=0, manufactured solution for rho1 and rho2 
    ! (to benchmark eqn1) Initializing rho1, rho2=Gaussian and rhototal has no-time dependence.
    !==================================================================================
 
    do j=lo(2),hi(2)
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1),hi(1)
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2
            w1  = alpha1*dexp(-rsq/(2.0d0*sigma**2))
            w2  =  delta*dexp(-beta*time)
            rhot = 1.0d0 + (molmass(2)*Dbar(3)/(molmass(1)*Dbar(1))-1.0d0)*w1
            rho(i,j,1) = rhot*w1
            rho(i,j,2) = rhot*w2 
            rho(i,j,3) = rhot-rho(i,j,1)-rho(i,j,2)
           
         end do
    end do

    case(6) 
    !=========================================================
    ! Test of thermodiffusion steady-state for 2 species 
    !=========================================================
 
    do j=lo(2),hi(2)
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1),hi(1)
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
   
            ! Solution to diff(c(y),y)=-K*c(y)*(1-c(y))
            ! Here K=grad(T)*S_T=0.15
            ! Height of domain H=32
            ! And average <rho1>=.4830852506
            rho(i,j,1) = 1.0d0/(1.0d0+0.1d0*exp(0.15d0*y))
            rho(i,j,2) = 1.0d0 - rho(i,j,1) 
   
         end do
      end do

   case default
      
      call bl_error("init_rho_2d: prob_type not supported")
      
    end select
   
  end subroutine init_rho_2d

  subroutine init_rho_3d(rho,ng_r,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng_r
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:) ! Last dimension for species 
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k,n
    real(kind=dp_t)  :: x,y,z,rsq,w1,w2,rhot,L(3),sum,rand,rho_loc,y1,r

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

    ! for specific box, now start loop over alloted cells
    select case (prob_type)
    
    case(0) 
    !================================================================================
    ! Thermodynamic equilibrium
    !================================================================================
 
    !$omp parallel do private(i,j,k)
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             rho(i,j,k,1:nspecies) = c_init(1,1:nspecies)

          end do
       end do
    end do
    !$omp end parallel do
    
    case(1) 
    !================================================================================
    ! Initializing rho's in concentric circle 
    !================================================================================
 
    !$omp parallel do private(i,j,k,x,y,z,rsq)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

               rsq = x**2 + y**2 + z**2
               if (rsq .lt. L(1)*L(2)*L(3)*0.001d0) then
                  rho(i,j,k,1:nspecies) = c_init(1,1:nspecies)
               else
                  rho(i,j,k,1:nspecies) = c_init(2,1:nspecies)
               end if
          
          end do
       end do
    end do
    !$omp end parallel do

    case(2) 
    !========================================================
    ! Initializing rho's with constant gradient
    !========================================================
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3),hi(3)
       z = prob_lo(3) + (dble(k)+half)*dx(3) 
       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2)
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1)

               rho(i,j,k,1:nspecies) = c_init(1,1:nspecies) + &
                   (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do
     end do
     !$omp end parallel do

     case(3) 
     !================================================================================
     ! Initializing rho's in Gaussian so as rhotot=constant=1.0. Here rho_exact = 
     ! e^(-r^2/4Dt)/(4piDt)^3/2, For norm, sigma/dx >2 (at t=0) & L/sigma < 8 (at t=t)
     !================================================================================
 
     !$omp parallel do private(i,j,k,x,y,z,rsq)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
        
              rsq = x**2 + y**2 + z**2
              rho(i,j,k,1) = dexp(-rsq/(4.0d0*Dbar(1)*time))/(4.0d0*M_PI*&
                             Dbar(1)*time)**1.5d0
              rho(i,j,k,2) = 1.0d0 - dexp(-rsq/(4.0d0*Dbar(1)*time))/(4.0d0*&
                             M_PI*Dbar(1)*time)**1.5d0
       
           end do
        end do
     end do
     !$omp end parallel do

     case(4)
     !==============================================================================
     ! Initializing rho1,rho2=Gaussian and rhot=space varying-constant 
     ! in time. Manufactured solution rho1_exact = e^(-r^2/4Dt-t/tau)/(4piDt)^(3/2)
     !==============================================================================
 
     !$omp parallel do private(i,j,k,x,y,z,rsq,rhot)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half) * dx(3) - half
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half) * dx(2) - half
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
              rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2 + (z-L(3)*half)**2
              rhot = 1.0d0 + alpha1*dexp(-rsq/(4.0d0*Dbar(1)))/(4.0d0*M_PI*&
                     Dbar(1))**1.5d0
           
              rho(i,j,k,1) = 1.0d0/(4.0d0*M_PI*Dbar(1)*time)**1.5d0*dexp(-rsq/&
                             (4.0d0*Dbar(1)*time) - time*beta)*rhot
              rho(i,j,k,2) = rhot - rho(i,j,k,1) 

           end do
        end do
     end do
     !$omp end parallel do

     case(5)
     !==================================================================================
     ! Initializing m2=m3, D12=D13 where Dbar(1)=D12, Dbar(2)=D13, Dbar(3)=D23, 
     ! Grad(w2)=0, manufactured solution for rho1 and rho2 
     !==================================================================================
 
     !$omp parallel do private(i,j,k,x,y,z,rsq,w1,w2,rhot)
     do k=lo(3),hi(3)
        z = prob_lo(3) + (dble(k)+half) * dx(3) - half
        do j=lo(2),hi(2)
           y = prob_lo(2) + (dble(j)+half) * dx(2) - half
           do i=lo(1),hi(1)
              x = prob_lo(1) + (dble(i)+half) * dx(1) - half

              rsq = (x-L(1)*half)**2 + (y-L(2)*half)**2 + (z-L(3)*half)**2
              w1  = alpha1*dexp(-rsq/(2.0d0*sigma**2))
              w2  =  delta*dexp(-beta*time)
              rhot = 1.0d0 + (molmass(2)*Dbar(3)/(molmass(1)*Dbar(1))-1.0d0)*w1
              rho(i,j,k,1) = rhot*w1
              rho(i,j,k,2) = rhot*w2
              rho(i,j,k,3) = rhot-rho(i,j,k,1)-rho(i,j,k,2)
           
              if(rho(i,j,k,1).lt.0.d0 .or. rho(i,j,k,2).lt.0.d0 .or. rho(i,j,k,3).lt.0.d0) then 
                 write(*,*), "rho1 / rho2 / rho3 is negative: STOP"
                 write(*,*), i, j, " w1=", w1, " w2=", w2, " rho1=",rho(i,j,k,1)," rho2=",&
                             rho(i,j,k,2), " rho3=",rho(i,j,k,3), " rhot=",rhot
              end if
 
          end do
       end do
    end do
    !$omp end parallel do

    case default
      
      call bl_error("init_rho_3d: prob_type not supported")
      
    end select
   
  end subroutine init_rho_3d

end module init_diffusion_module
