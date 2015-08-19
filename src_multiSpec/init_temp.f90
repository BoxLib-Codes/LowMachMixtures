module init_temp_module

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
  use probin_common_module, only: prob_lo, prob_hi
  use probin_multispecies_module, only: T_init, nspecies, sigma, temp_type
 
  implicit none

  private

  public :: init_Temp
  
contains

  subroutine init_Temp(Temp,dx,time,the_bc_level)

    type(multifab) , intent(inout) :: Temp(:)            
    real(kind=dp_t), intent(in   ) :: dx(:,:)           
    real(kind=dp_t), intent(in   ) :: time 
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(Temp(1)%dim), hi(Temp(1)%dim)
    integer                        :: dm, ng, i, n, nlevs
    real(kind=dp_t), pointer       :: dp1(:,:,:,:)  ! pointer for Temp 

    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_Temp")

    dm = Temp(1)%dim
    ng = Temp(1)%ng
    nlevs = size(Temp,1)

    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(Temp(n))
          dp1 => dataptr(Temp(n),i)
          lo  = lwb(get_box(Temp(n),i))
          hi  = upb(get_box(Temp(n),i))
          
          select case(dm)
          case (2)
             call init_Temp_2d(dp1(:,:,1,1),ng,lo,hi,dx(n,:),time)
          case (3)
             call init_Temp_3d(dp1(:,:,:,1),ng,lo,hi,dx(n,:),time)
          end select
       end do

       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(Temp(n))

       ! fill non-periodic domain boundary ghost cells
       call multifab_coefbc(Temp(n),1,1,the_bc_level(n))

    end do

    call destroy(bpt)

  end subroutine init_Temp

  subroutine init_Temp_2d(Temp,ng,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:)  
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,rsq,L(2)
 
    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    
    ! for specific box, now start loop over alloted cells     
    select case (temp_type)

    case(0) 
    !=========================================================================
    ! Thermodynamic equilibrium
    !=========================================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
 
            Temp(i,j) = T_init(1)
            ! These were used for testing thermodiffusion only
            if(.false.) Temp(i,j) = 1.0d0 + 0.01d0*cos(2.0d0*M_PI*x/L(1))*sin(2.0d0*M_PI*y/L(1))
            if(.false.) Temp(i,j) = 1.0d0 + 0.01d0*sin(2.0d0*M_PI*y/L(1))

       end do
    end do  
    
    case(1) 
    !=========================================================================
    ! Initializing T in concentric circle at (Lx/2,Ly/2) with radius^2=0.1*L(1)*L(2)
    !=========================================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
      
          ! temperature distribution follows the density 
          rsq = x**2 + y**2
          if (rsq .lt. L(1)*L(2)*0.1d0) then
              Temp(i,j) = T_init(1)
          else
              Temp(i,j) = T_init(2)
          end if
    
        end do
    end do
  
    case(2,6) 
    !========================================================
    ! Initializing T with constant gradient along y axes
    !========================================================
    do j=lo(2)-ng,hi(2)+ng
       y = prob_lo(2) + (dble(j)+half)*dx(2) 
       do i=lo(1)-ng,hi(1)+ng
          x = prob_lo(1) + (dble(i)+half)*dx(1) 
      
          ! linear gradient in y direction
          Temp(i,j) = T_init(1) + (T_init(2) - T_init(1))*(y-prob_lo(2))/L(2)

       end do
    end do

    case default

    do j=lo(2)-ng,hi(2)+ng
         y = prob_lo(2) + (dble(j)+half) * dx(2) - half
         do i=lo(1)-ng,hi(1)+ng
            x = prob_lo(1) + (dble(i)+half) * dx(1) - half
        
            Temp(i,j)  = T_init(1) 

         end do
    end do

   end select
   
  end subroutine init_Temp_2d

  subroutine init_Temp_3d(Temp,ng,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng
    real(kind=dp_t)  :: Temp(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:)  
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k
    real(kind=dp_t)  :: x,y,z,rsq,L(3)

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length

    ! for specific box, now start loop over alloted cells     
    select case (temp_type)
    
    case(0) 
    !================================================================================
    ! Thermodynamic equilibrium
    !================================================================================
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng
             
             Temp(i,j,k) = T_init(1)

          end do
       end do
    end do
    !$omp end parallel do
    
    case(1) 
    !================================================================================
    ! Initializing temperature in concentric circle 
    !================================================================================
  
    !$omp parallel do private(i,j,k,x,y,z,rsq)
    do k=lo(3)-ng,hi(3)+ng
       z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
             
             !temperature distribution follows the density 
             rsq = x**2 + y**2 + z**2
             if (rsq .lt. L(1)*L(2)*L(3)*0.001d0) then
                Temp(i,j,k) = T_init(1)
             else
                Temp(i,j,k) = T_init(2)
             end if
          
          end do
       end do
    end do
    !$omp end parallel do

    case(2,6) 
    !========================================================
    ! Initializing T with constant gradient along y direction
    !========================================================
 
    !$omp parallel do private(i,j,k,x,y,z)
    do k=lo(3)-ng,hi(3)+ng
       z = prob_lo(3) + (dble(k)+half)*dx(3) 
       do j=lo(2)-ng,hi(2)+ng
          y = prob_lo(2) + (dble(j)+half)*dx(2) 
          do i=lo(1)-ng,hi(1)+ng
             x = prob_lo(1) + (dble(i)+half)*dx(1) 

             ! linear gradient in y direction for temperature 
             Temp(i,j,k) = T_init(1) + (T_init(2) - T_init(1))*(y-prob_lo(2))/L(2)
 
          end do
       end do
    end do
    !$omp end parallel do

    case default

    !$omp parallel do private(i,j,k)
    do k=lo(3)-ng,hi(3)+ng
       do j=lo(2)-ng,hi(2)+ng
          do i=lo(1)-ng,hi(1)+ng

             Temp(i,j,k) = T_init(1)

          end do
       end do
    end do
    !$omp end parallel do

   end select
   
  end subroutine init_Temp_3d

end module init_temp_module
