module external_force_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use multifab_fill_ghost_module
  use probin_multispecies_module, only: Dbar, sigma, beta, alpha1, delta
  use probin_common_module, only: prob_lo, prob_hi, prob_type, molmass

  implicit none

  private

  public :: external_source

contains

  subroutine external_source(mla,rho,fluxdiv,dx,time)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(inout) :: fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),time

    ! local variables
    integer :: i,n,ng,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointer for rho,fluxdiv
    real(kind=dp_t), pointer :: dp(:,:,:,:)
    real(kind=dp_t), pointer :: dp1(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"external_source")

    dm    = mla%dim         ! dimensionality
    ng    = fluxdiv(1)%ng   ! number of ghost cells
    nlevs = mla%nlevel      ! number of levels

    do n=1,nlevs
       do i = 1, nfabs(rho(n))
          dp  => dataptr(rho(n),i)
          dp1 => dataptr(fluxdiv(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call external_source_2d(dp(:,:,1,:),dp1(:,:,1,:),&
                                     ng,lo,hi,dx(n,:),time)
          case (3)
             call external_source_3d(dp(:,:,:,:),dp1(:,:,:,:),&
                                     ng,lo,hi,dx(n,:),time)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine external_source

  subroutine external_source_2d(rho,fluxdiv,ng,lo,hi,dx,time)

    integer          :: lo(:),hi(:),ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t)  :: fluxdiv(lo(1)-ng:,lo(2)-ng:,:)
    real(kind=dp_t)  :: dx(:),time

    ! local variables
    integer          :: i,j
    real(kind=dp_t)  :: x,y,r,L(2),D12,D23,r_temp,r_temp1

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length
    D12    = Dbar(1)
    D23    = Dbar(3)
   
    select case (abs(prob_type))
    
    !==== for m1 = m2 = m and also for m1 != m2 != m, 2 species ====!
    case(4) 

     ! for specific box, now start loops over alloted cells    
     do j=lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
        do i=lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
 
             r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2)
               
             r_temp = -(dexp(-(beta*time) - (r**2*(1.0d0 + time))/(4.0d0*D12*time))*(16.0d0*beta*D12**2*&
                      dexp(r**2/(4.0d0*D12))*M_PI*time + alpha1*(r**2 + 4.0d0*beta*D12*time)))/(64.0d0*&
                      D12**3*M_PI**2*time**2)
  
             fluxdiv(i,j,1) = fluxdiv(i,j,1) + r_temp 
             fluxdiv(i,j,2) = fluxdiv(i,j,2) - r_temp

        end do
     end do
  

    !==== for m2=m3, D12=D13 and Grad(w2)=0, 3 species ====! 
    case(5) 
     
     ! for specific box, now start loops over alloted cells    
     do j=lo(2), hi(2)
        y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
        do i=lo(1), hi(1)
           x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0
 
             r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2)
             r_temp  = (alpha1*(2.0d0*alpha1*(D12*molmass(1) - D23*molmass(2))*(r - sigma)*(r +& 
                       sigma) - D12*dexp(r**2/(2.0d0*sigma**2))*molmass(1)*(r**2 - 2.0d0*sigma**2)))/&
                       (dexp(r**2/sigma**2)*molmass(1)*sigma**4)

             r_temp1 = (dexp(-r**2/(2.0d0*sigma**2) - beta*time)*(-(beta*D12*delta*dexp(r**2/(2.0d0*&
                       sigma**2))*molmass(1)*sigma**4) + alpha1*delta*(D12*(D12 - D23)*molmass(1)*&
                       r**2 + 2.0d0*D12*(-D12 + D23)*molmass(1)*sigma**2 + beta*(D12*molmass(1) -& 
                       D23*molmass(2))*sigma**4)))/(D12*molmass(1)*sigma**4)
            
             fluxdiv(i,j,1) = fluxdiv(i,j,1) + r_temp 
             fluxdiv(i,j,2) = fluxdiv(i,j,2) + r_temp1
             fluxdiv(i,j,3) = fluxdiv(i,j,3) - (r_temp + r_temp1)
              
             !== to benchmark eqn1 ==!
             !r_temp = (alpha1*D12*(-(dexp(r**2/(2.0d0*sigma**2))*(r**2 - 2.0d0*sigma**2)) +& 
             !         2.0d0*beta*(-r**2 + sigma**2)))/(dexp(r**2/sigma**2)*sigma**4)  
             !fluxdiv(i,j,1) = fluxdiv(i,j,1) + r_temp 
             !fluxdiv(i,j,2) = fluxdiv(i,j,2) - r_temp
             !=======================!

        end do
     end do
 
    end select 

  end subroutine external_source_2d

  subroutine external_source_3d(rho,fluxdiv,ng,lo,hi,dx,time)

    integer          :: lo(:),hi(:),ng
    real(kind=dp_t)  :: rho(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t)  :: fluxdiv(lo(1)-ng:,lo(2)-ng:,lo(3)-ng:,:)
    real(kind=dp_t)  :: dx(:),time

    ! local variables
    integer          :: i,j,k
    real(kind=dp_t)  :: x,y,z,r,L(3),D12,D23,r_temp,r_temp1

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    D12    = Dbar(1)
    D23    = Dbar(3)
    
    select case (abs(prob_type))
    
    !==== for m1 = m2 = m and m1 != m2 != m, 2 species ====!
    case(4)
     
     ! for specific box, now start loops over alloted cells    
     do k=lo(3), hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
        do j=lo(2), hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
           do i=lo(1), hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0           
   
              r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2)
              
              r_temp = (dexp(-(beta*time) - (r**2*(1.0d0 + time))/(4.0d0*D12*time))*&
                       (0.02244839026564582d0*beta*D12**3.0d0*dexp(r**2/(4.0d0*D12))*(D12*&
                       time)**3.5d0 + alpha1*D12**1.5d0*(0.0001259825563796855d0*r**2*(D12*time)**2.5d0 +& 
                       0.000503930225518742d0*beta*(D12*time)**3.5d0)))/(D12**3.0d0*(D12*time)**5.0d0) 
 
              fluxdiv(i,j,k,1) = fluxdiv(i,j,k,1) - r_temp 
              fluxdiv(i,j,k,2) = fluxdiv(i,j,k,2) + r_temp 
 
           end do
        end do
     end do

    !==== for m2=m3, D12=D13 and Grad(w2)=0, 3 species ====! 
    case(5) 
    
     do k=lo(3), hi(3)
        z = prob_lo(3) + (dble(k)+0.5d0) * dx(3) - 0.5d0
        do j=lo(2), hi(2)
           y = prob_lo(2) + (dble(j)+0.5d0) * dx(2) - 0.5d0
           do i=lo(1), hi(1)
              x = prob_lo(1) + (dble(i)+0.5d0) * dx(1) - 0.5d0           
   
              r = sqrt((x-L(1)*0.5d0)**2 + (y-L(2)*0.5d0)**2 + (z-L(3)*0.5d0)**2)
              
               r_temp = (alpha1*(-(D12*dexp(r**2/(2.0d0*sigma**2))*molmass(1)*(r**2 - 3.0d0*sigma**2)) +& 
                        alpha1*(D12*molmass(1) - D23*molmass(2))*(2.0d0*r**2 - 3.0d0*sigma**2)))/&
                        (dexp(r**2/sigma**2)*molmass(1)*sigma**4) 
                       
              r_temp1 = (dexp(-r**2/(2.0d0*sigma**2) - beta*time)*(-(beta*D12*delta*dexp(r**2/(2.0d0*&
                        sigma**2))*molmass(1)*sigma**4) + alpha1*delta*(D12*(D12 - D23)*molmass(1)*&
                        r**2 + 3.0d0*D12*(-D12 + D23)*molmass(1)*sigma**2 + beta*(D12*molmass(1) -& 
                        D23*molmass(2))*sigma**4)))/(D12*molmass(1)*sigma**4)       
 
              fluxdiv(i,j,k,1) = fluxdiv(i,j,k,1) + r_temp 
              fluxdiv(i,j,k,2) = fluxdiv(i,j,k,2) + r_temp1 
              fluxdiv(i,j,k,3) = fluxdiv(i,j,k,3) - (r_temp + r_temp1) 

           end do
        end do
     end do


    end select 
 
  end subroutine external_source_3d

end module external_force_module 
