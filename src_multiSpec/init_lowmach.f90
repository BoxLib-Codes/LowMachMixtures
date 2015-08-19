module init_lowmach_module

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
  use probin_common_module, only: prob_lo, prob_hi, prob_type, k_B, grav, &
                                  molmass, rhobar, smoothing_width, u_init, n_cells
  use probin_multispecies_module, only: alpha1, beta, delta, sigma, Dbar, Dtherm, &
                                        c_init, nspecies, T_init
 
  implicit none

  private

  public :: init_rho_and_umac  ! used in low Mach code; initialize c first then convert to rho

  ! IMPORTANT: In the diffusion only code (init_rho), c_init specifies initial values for DENSITY
  ! In the low-Mach code (init_rho_and_umac), c_init specifies initial MASS FRACTIONS
  ! (should sum to unity!... but we overwrite the final concentration so sum(c_i)=1 before computing rho)
  ! The density follows from the EOS in the LM case so it cannot be specified
  ! Same applies to boundary conditions

  ! prob_types codes for init_lowmach:

  !=============================================================
  ! case 1:
  ! bubble with radius = 1/4 of domain in x
  ! c=c_init(1,:) inside, c=c_init(2,:) outside
  ! can be discontinous or smooth depending on smoothing_width

  !=========================================================
  ! case 2:
  ! constant concentration gradient along y
  ! c=c_init(1,:) on bottom, c=c_init(2,:) on top

  !=========================================================
  ! case 3:
  ! 1 fluid on top of another
  ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
  ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
  !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
  ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
  ! smoothing width of -2 is a sinusoidal perturbation
  ! x-vel = u_init(1) below centerline, u_init(2) above centerline

  !=========================================================
  ! case 4:
  ! not defined
  ! external_source has an analytic solution, not sure what it was for
  ! I think it's a remnant from test_diffusion

  !=========================================================
  ! case 5:
  ! not defined
  ! external_source has an analytic solution, not sure what it was for
  ! I think it's a remnant from test_diffusion

  !=========================================================
  ! case 6:
  ! Two Gaussian bubbles with different centers
  ! c(1) peak is c_init(1,1) 
  ! c(2) peak is c_init(2,2) 

  !=========================================================
  ! case 7:
  ! not defined

  !=========================================================
  ! case 8:
  ! not defined

  !=========================================================
  ! case 9:
  ! not defined
  ! mixture_properties_mass_local computes Dbar's using water/glycerol
  ! compute_eta uses water/glycerol

  !=========================================================
  ! case 10:
  ! not defined

  !=========================================================
  ! case 11:
  ! Discontinuous square in the central 25% of domain
  ! c=c_init(1,:) inside; c=c_init(2,:) outside

  !=========================================================
  ! case 12:
  ! Gaussian bubble centered in domain
  ! c=c_init(1,:) inside; c=c_init(2,:) outside
  ! lo- and hi-y walls move with prescribed velocity,
  ! see inhomogeneous_bc_val.f90
  ! compute_eta uses linear profile in rho if prob_type < 0

  !=========================================================
  ! case 13:
  ! stratified multispecies due to barodiffusion
  ! approximate analytical steady solution
  ! assumes the final species is the dominant component

  !=========================================================
  ! case 14:
  ! stratified multispecies due to thermodiffusion
  ! approximate analytical steady solution
  ! assumes nspecies=3
  ! assumes the final species is the dominant component


contains

  subroutine init_rho_and_umac(mla,rho,umac,dx,time,the_bc_level)

    ! initialize rho_i and umac in the valid region
    ! we first initialize c_i in the valid region
    ! then enforce that sum(c_i)=1 by overwriting the final concentration,
    ! and then use the EOS to compute rho_i

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    real(kind=dp_t), intent(in   ) :: time
    type(bc_level) , intent(in   ) :: the_bc_level(:)
 
    ! local variables
    integer                        :: lo(mla%dim), hi(mla%dim)
    integer                        :: dm, ng_c, ng_u, ng_r, i, n, nlevs
    real(kind=dp_t), pointer       :: dp(:,:,:,:)   ! pointer for rho (last dim:nspecies)   
    real(kind=dp_t), pointer       :: up(:,:,:,:)   ! pointers for mac velocities
    real(kind=dp_t), pointer       :: vp(:,:,:,:)
    real(kind=dp_t), pointer       :: wp(:,:,:,:)
    real(kind=dp_t), pointer       :: rp(:,:,:,:)

    type(multifab) :: conc(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "init_rho_and_umac")

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       call multifab_build(conc(n),mla%la(n),nspecies,rho(n)%ng)
    end do

    ng_u = umac(1,1)%ng
    ng_c = conc(1)%ng
    ng_r = rho(1)%ng

    ! assign values of parameters for the gaussian rho, rhototal
    alpha1 = 0.5d0 
    beta   = 0.1d0 
    delta  = 0.5d0 
    sigma  = (prob_hi(1)-prob_lo(1))/10.0d0  ! variance of gaussian distribution

    ! looping over boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(conc(n),i)
          rp => dataptr(rho(n),i)
          up => dataptr(umac(n,1),i)
          vp => dataptr(umac(n,2),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case (dm)
          case (2)
             call init_rho_and_umac_2d(dp(:,:,1,:),ng_c,rp(:,:,1,:),ng_r, &
                                       up(:,:,1,1),vp(:,:,1,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          case (3)
             wp => dataptr(umac(n,3),i)
             call init_rho_and_umac_3d(dp(:,:,:,:),ng_c,rp(:,:,:,:),ng_r, &
                                       up(:,:,:,1),vp(:,:,:,1),wp(:,:,:,1),ng_u, &
                                       lo,hi,dx(n,:),time)
          end select
       end do
    end do

    do n=1,nlevs
       call multifab_destroy(conc(n))
    end do

    call destroy(bpt)

  end subroutine init_rho_and_umac

  subroutine init_rho_and_umac_2d(c,ng_c,rho,ng_r,u,v,ng_u,lo,hi,dx,time)

    integer          :: lo(2), hi(2), ng_c, ng_u, ng_r
    real(kind=dp_t)  ::   c(lo(1)-ng_c:,lo(2)-ng_c:,:)
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local varables
    integer          :: i,j,n,seed(12)
    real(kind=dp_t)  :: x,y,rad,L(2),sum,r,r1,r2,y1,c_loc
    real(kind=dp_t)  :: gradToverT,m_e
 
    real(kind=dp_t)  :: random

    L(1:2) = prob_hi(1:2)-prob_lo(1:2) ! Domain length

    seed = n_cells(1)*lo(2) + lo(1) + 1
    call random_seed(put=seed(1:12))

    select case (abs(prob_type))
    
    case (1)

       !=============================================================
       ! bubble with radius = 1/4 of domain in x
       ! c=c_init(1,:) inside, c=c_init(2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=============================================================
 
       u = 0.d0
       v = 0.d0

       rad = L(1)/4.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))
       
             r = sqrt(x**2 + y**2)

             if (smoothing_width .eq. 0) then

                ! discontinuous interface
                if (r .lt. rad) then
                   c(i,j,1:nspecies) = c_init(1,1:nspecies)
                else
                   c(i,j,1:nspecies) = c_init(2,1:nspecies)
                end if

             else

                ! smooth interface
                c(i,j,1:nspecies-1) = c_init(1,1:nspecies-1) + &
                     (c_init(2,1:nspecies-1) - c_init(1,1:nspecies-1))* &
                     0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

             end if
    
          end do
       end do

    case (2) 

       !=========================================================
       ! constant concentration gradient along y
       ! c=c_init(1,:) on bottom, c=c_init(2,:) on top
       !=========================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + (dble(j)+half)*dx(2) 
          do i=lo(1),hi(1)
             x = prob_lo(1) + (dble(i)+half)*dx(1) 

             ! linear gradient in mass fractions
             c(i,j,1:nspecies) = c_init(1,1:nspecies) + & 
                  (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

          end do
       end do

    case (3)

       !=============================================================
       ! 1 fluid on top of another
       ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
       ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
       !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
       ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
       ! smoothing width of -2 is a sinusoidal perturbation
       ! x-vel = u_init(1) below centerline, u_init(2) above centerline
       !=============================================================

       u = 0.d0
       v = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       ! x-velocity = u_init(1) below centerline
       !              u_init(2) above centerline
       do j=lo(2),hi(2)
          y = prob_lo(2) + (j+0.5d0)*dx(2)
          if (y .lt. y1) then
             u(:,j) = u_init(1)
          else
             u(:,j) = u_init(2)
          end if
       end do

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then
          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   call random_number(random)
                   c_loc = abs(smoothing_width)*random
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do n=1,nspecies
                c_loc = c_init(1,n) + (c_init(2,n)-c_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,n) = c_loc
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do i=lo(1),hi(1)
                   x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                   c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                   do n=1,nspecies
                      c(i,j,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_2d: smoothing_width not compatible with prob_type")

       end if

    case (6)

       !=============================================================
       ! Two Gaussian bubbles with different centers
       ! c(1) peak is c_init(1,1) 
       ! c(2) peak is c_init(2,2) 
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)

             r1 = sqrt ((x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 + (y-(0.6d0*prob_lo(2)+0.4d0*prob_hi(2)))**2)
             r2 = sqrt ((x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 + (y-(0.4d0*prob_lo(2)+0.6d0*prob_hi(2)))**2)
             

             ! set c using Gaussian bump
             c(i,j,1) = c_init(1,1)*exp(-75.d0*r1**2)
             c(i,j,2) = c_init(2,2)*exp(-75.d0*r2**2)

          enddo
       enddo

    case (11)

       !=============================================================
       ! Discontinuous square in the central 25% of domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             ! initialize c to a square region
             if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                 j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1) then
                c(i,j,1:nspecies) = c_init(1,1:nspecies)
             else
                c(i,j,1:nspecies) = c_init(2,1:nspecies)
             end if

          enddo
       enddo

    case (12)

       !=============================================================
       ! Gaussian bubble centered in domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       ! lo- and hi-y walls move with prescribed velocity,
       ! see inhomogeneous_bc_val.f90
       !=============================================================

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
          do i=lo(1),hi(1)
             x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

             r = sqrt (x**2 + y**2)

             ! set c using Gaussian bump
             c(i,j,1:nspecies-1) = c_init(1,1:nspecies-1)*exp(-75.d0*r**2)

          enddo
       enddo

    case (13)

       !=============================================================
       ! stratified multispecies due to barodiffusion
       ! approximate analytical steady solution
       ! assumes the final species is the dominant component
       !=============================================================

       u = 0.d0
       v = 0.d0

       do n=1,nspecies-1
          m_e = (rhobar(nspecies)/rhobar(n) - 1.d0)*molmass(n)

          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
                c(i,j,n) = c_init(1,n)*exp(-m_e*grav(2)*y/(k_B*T_init(1)))
             enddo
          enddo
       enddo

    case (14)

       !=============================================================
       ! stratified multispecies due to thermodiffusion
       ! approximate analytical steady solution
       ! assumes nspecies=3
       ! assumes the final species is the dominant component
       !=============================================================

       if (nspecies .ne. 3) then
          call bl_error("prob_type=14 requires nspecies=3")
       end if

       gradToverT = (T_init(2)-T_init(1))/(T_init(1)*(prob_hi(2)-prob_lo(2)))

       u = 0.d0
       v = 0.d0

       do j=lo(2),hi(2)
          y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
          do i=lo(1),hi(1)
             
             c(i,j,1) = c_init(1,1)*exp((Dtherm(1)-Dtherm(3))*gradToverT*y/Dbar(2))
             c(i,j,2) = c_init(1,2)*exp((Dtherm(2)-Dtherm(3))*gradToverT*y/Dbar(3))

          enddo
       enddo

    case default

       call bl_error("Desired prob_type not supported in 2D")

    end select

    ! set final c_i such that sum(c_i) = 1 to within roundoff
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0
          do n=1,nspecies-1
             sum = sum + c(i,j,n)
          end do
          c(i,j,nspecies) = 1.d0 - sum

       end do
    end do

    ! compute rho using the eos
    do j=lo(2),hi(2)
       do i=lo(1),hi(1)

          sum = 0.d0
          do n=1,nspecies
             ! sum represents rhoinv
             sum = sum + c(i,j,n)/rhobar(n)
          end do
          rho(i,j,1:nspecies) = c(i,j,1:nspecies)/sum

       end do
    end do

  end subroutine init_rho_and_umac_2d

  subroutine init_rho_and_umac_3d(c,ng_c,rho,ng_r,u,v,w,ng_u,lo,hi,dx,time)
    
    integer          :: lo(3), hi(3), ng_c, ng_u, ng_r
    real(kind=dp_t)  ::   c(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
    real(kind=dp_t)  :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)
    real(kind=dp_t)  ::   u(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   v(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  ::   w(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t)  :: dx(:)
    real(kind=dp_t)  :: time 
 
    ! local variables
    integer          :: i,j,k,n,seed(12)
    real(kind=dp_t)  :: x,y,z,rad,L(3),sum,c_loc,y1,r,r1,r2,m_e,gradToverT

    real(kind=dp_t) :: random

    L(1:3) = prob_hi(1:3)-prob_lo(1:3) ! Domain length
    
    seed = n_cells(1)*n_cells(2)*lo(3) + n_cells(1)*lo(2) + lo(1)
    call random_seed(put=seed(1:12))

    select case (abs(prob_type))

    case (1) 

       !=============================================================
       ! bubble with radius = 1/4 of domain in x
       ! c=c_init(1,:) inside, c=c_init(2,:) outside
       ! can be discontinous or smooth depending on smoothing_width
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       rad = L(1)/4.d0

       !$omp parallel do private(i,j,k,x,y,z,r)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+half)*dx(3) - half*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+half)*dx(2) - half*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+half)*dx(1) - half*(prob_lo(1)+prob_hi(1))

                r = sqrt(x**2 + y**2 + z**2)

                if (smoothing_width .eq. 0) then

                   ! discontinuous interface
                   if (r .lt. rad) then
                      c(i,j,k,1:nspecies) = c_init(1,1:nspecies)
                   else
                      c(i,j,k,1:nspecies) = c_init(2,1:nspecies)
                   end if

                else

                   ! smooth interface
                   c(i,j,k,1:nspecies-1) = c_init(1,1:nspecies-1) + &
                        (c_init(2,1:nspecies-1) - c_init(1,1:nspecies-1))* &
                        0.5d0*(1.d0 + tanh((r-rad)/(smoothing_width*dx(1))))

                end if

             end do
          end do
       end do
       !$omp end parallel do

    case (2) 

       !=========================================================
       ! constant concentration gradient along y
       ! c=c_init(1,:) on bottom, c=c_init(2,:) on top
       !=========================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       !$omp parallel do private(i,j,k,x,y,z)
       do k=lo(3),hi(3)
          z = prob_lo(3) + (dble(k)+half)*dx(3) 
          do j=lo(2),hi(2)
             y = prob_lo(2) + (dble(j)+half)*dx(2)
             do i=lo(1),hi(1)
                x = prob_lo(1) + (dble(i)+half)*dx(1)

                c(i,j,k,1:nspecies) = c_init(1,1:nspecies) + &
                     (c_init(2,1:nspecies) - c_init(1,1:nspecies))*(y-prob_lo(2))/L(2)

             end do
          end do
       end do
       !$omp end parallel do

    case (3) 

       !=============================================================
       ! 1 fluid on top of another
       ! c = c_init(1,:) on bottom; c = c_init(2,:) on top
       ! smoothing_width > 0 is a tanh smoothed interface where smoothing width is approx the # of grid 
       !   cells and then c = rand*c_init(1,:) + (1-rand)*c_init(2,:)
       ! smoothing_width between 0 and -1 is random perturbation where rand = abs(smoothing_width)*rand()
       ! smoothing width of -2 is a sinusoidal perturbation
       ! x-vel = u_init(1) below centerline, u_init(2) above centerline
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       ! middle of domain
       y1 = (prob_lo(2)+prob_hi(2)) / 2.d0

       ! x-velocity = u_init(1) below centerline
       !              u_init(2) above centerline
       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                u(:,j,k) = u_init(1)
             else
                u(:,j,k) = u_init(2)
             end if
          end do
       end do

       if (smoothing_width .le. 0.d0 .and. smoothing_width .ge. -1.d0) then

          ! discontinuous version with random perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(2,n)
                end do
             end if

             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      call random_number(random)
                      c_loc = abs(smoothing_width)*random
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                      end do
                   end do
                end do
             end if

          end do

       else if (smoothing_width .gt. 0.d0) then

          ! smoothed version
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0) - y1
             do n=1,nspecies
                c_loc = c_init(1,n) + (c_init(2,n)-c_init(1,n))*0.5d0*(tanh(y/(smoothing_width*dx(2)))+1.d0)
                c(lo(1):hi(1),j,lo(3):hi(3),n) = c_loc
             end do
          end do

       else if (smoothing_width .eq. -2.d0) then

          ! discontinuous version with sinusoidal perturbation
          do j=lo(2),hi(2)
             y = prob_lo(2) + (j+0.5d0)*dx(2)
             if (y .lt. y1) then
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(1,n)
                end do
             else
                do n=1,nspecies
                   c(lo(1):hi(1),j,lo(3):hi(3),n) = c_init(2,n)
                end do
             end if
             if (j .eq. n_cells(2)/2) then
                do k=lo(3),hi(3)
                   do i=lo(1),hi(1)
                      x = prob_lo(1) + (dble(i)+0.5d0)*dx(1)
                      c_loc = 0.5d0*(cos(4.d0*M_PI*x/L(1))+1.d0)
                      do n=1,nspecies
                         c(i,j,k,n) = c_loc*(c_init(1,n)) + (1.d0-c_loc)*c_init(2,n)
                      end do
                   end do
                end do
             end if

          end do

       else

          call bl_error("init_rho_and_umac_3d: smoothing_width not compatible with prob_type")

       end if

    case (6)

       !=============================================================
       ! Two Gaussian bubbles with different centers
       ! c(1) peak is c_init(1,1)
       ! c(2) peak is c_init(2,2)
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0
       
       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0)
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0)
                
                r1 = sqrt (   (x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 &
                            + (y-(0.6d0*prob_lo(2)+0.4d0*prob_hi(2)))**2 &
                            + (z-(0.5d0*prob_lo(3)+0.5d0*prob_hi(3)))**2)
                r2 = sqrt (   (x-(0.5d0*prob_lo(1)+0.5d0*prob_hi(1)))**2 &
                            + (y-(0.4d0*prob_lo(2)+0.6d0*prob_hi(2)))**2 &
                            + (z-(0.5d0*prob_lo(3)+0.5d0*prob_hi(3)))**2)
                
                ! set c using Gaussian bump
                c(i,j,k,1) = c_init(1,1)*exp(-75.d0*r1**2)
                c(i,j,k,2) = c_init(2,2)*exp(-75.d0*r2**2)
                
             enddo
          enddo
       enddo

    case (11)

       !=============================================================
       ! Discontinuous square in the central 25% of domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)

                ! initialize c to a square region
                if (i .ge. n_cells(1)/4 .and. i .le. 3*n_cells(1)/4-1 .and. &
                    j .ge. n_cells(2)/4 .and. j .le. 3*n_cells(2)/4-1 .and. &
                    k .ge. n_cells(3)/4 .and. k .le. 3*n_cells(3)/4-1) then
                   c(i,j,k,1:nspecies) = c_init(1,1:nspecies)
                else
                   c(i,j,k,1:nspecies) = c_init(2,1:nspecies)
                end if

             enddo
          enddo
       enddo

    case (12)

       !=============================================================
       ! Gaussian bubble centered in domain
       ! c=c_init(1,:) inside; c=c_init(2,:) outside
       ! lo- and hi-y walls move with prescribed velocity,
       ! see inhomogeneous_bc_val.f90
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          z = prob_lo(3) + dx(3) * (dble(k)+0.5d0) - 0.5d0*(prob_lo(3)+prob_hi(3))
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2) * (dble(j)+0.5d0) - 0.5d0*(prob_lo(2)+prob_hi(2))
             do i=lo(1),hi(1)
                x = prob_lo(1) + dx(1) * (dble(i)+0.5d0) - 0.5d0*(prob_lo(1)+prob_hi(1))

                r = sqrt (x**2 + y**2 + z**2)

                ! set c using Gaussian bump
                c(i,j,k,1:nspecies-1) = c_init(1,1:nspecies-1)*exp(-75.d0*r**2)

             enddo
          enddo
       enddo

    case (13)

       !=============================================================
       ! stratified multispecies due to barodiffusion
       ! approximate analytical steady solution
       ! assumes the final species is the dominant component
       !=============================================================

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do n=1,nspecies-1
          m_e = (rhobar(nspecies)/rhobar(n) - 1.d0)*molmass(n)
          do k=lo(3),hi(3)
             do j=lo(2),hi(2)
                y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
                do i=lo(1),hi(1)
                   c(i,j,k,n) = c_init(1,n)*exp(-m_e*grav(2)*y/(k_B*T_init(1)))
                enddo
             enddo
          enddo
       enddo

    case (14)

       !=============================================================
       ! stratified multispecies due to thermodiffusion
       ! approximate analytical steady solution
       ! assumes nspecies=3
       ! assumes the final species is the dominant component
       !=============================================================

       if (nspecies .ne. 3) then
          call bl_error("prob_type=14 requires nspecies=3")
       end if

       gradToverT = (T_init(2)-T_init(1))/(T_init(1)*(prob_hi(2)-prob_lo(2)))

       u = 0.d0
       v = 0.d0
       w = 0.d0

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             y = prob_lo(2) + dx(2)*(dble(j)+0.5d0)
             do i=lo(1),hi(1)
             
                c(i,j,k,1) = c_init(1,1)*exp((Dtherm(1)-Dtherm(3))*gradToverT*y/Dbar(2))
                c(i,j,k,2) = c_init(1,2)*exp((Dtherm(2)-Dtherm(3))*gradToverT*y/Dbar(3))

             enddo
          enddo
       enddo

    case default

       call bl_error("Desired prob_type not supported in 3D")

    end select

    ! set final c_i such that sum(c_i) = 1
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             sum = 0
             do n=1,nspecies-1
                sum = sum + c(i,j,k,n)
             end do
             c(i,j,k,nspecies) = 1.d0 - sum

          end do
       end do
    end do

    ! compute rho using the eos
    do k=lo(3),hi(3)
       do j=lo(2),hi(2)
          do i=lo(1),hi(1)

             sum = 0.d0
             do n=1,nspecies
                ! sum represents rhoinv
                sum = sum + c(i,j,k,n)/rhobar(n)
             end do
             rho(i,j,k,1:nspecies) = c(i,j,k,1:nspecies)/sum

          end do
       end do
    end do

  end subroutine init_rho_and_umac_3d

end module init_lowmach_module
