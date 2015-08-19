module stochastic_m_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use BoxLibRNGs
  use multifab_physbc_stag_module
  use multifab_fill_random_module
  use multifab_filter_module
  use sum_momenta_module
  use convert_m_to_umac_module
  use convert_stag_module
  use probin_common_module , only: visc_coef, variance_coef_mom, k_B, &
                                   stoch_stress_form, filtering_width

  implicit none

  interface add_m_fluctuations
     module procedure add_m_fluctuations_1
     module procedure add_m_fluctuations_2
  end interface

  private

  public :: stochastic_m_fluxdiv, fill_m_stochastic, &
       init_m_stochastic, destroy_m_stochastic, add_m_fluctuations

  ! Stochastic fluxes for momentum are generated on:
  ! -cell-centered grid for diagonal components
  ! -node-centered (2D) or edge-centered (3D) grid for off-diagonal components
  type(multifab), allocatable, save :: mflux_cc(:,:), mflux_nd(:,:), mflux_ed(:,:,:)

  integer, save :: n_rngs ! how many random number stages

contains

  ! Note that here we *increment* stoch_m_force so it must be initialized externally!
  subroutine stochastic_m_fluxdiv(mla,the_bc_level,stoch_m_force,eta,eta_ed, &
                                  temperature,temperature_ed,dx,dt,weights)
    
    type(ml_layout), intent(in   ) :: mla
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    type(multifab) , intent(inout) :: stoch_m_force(:,:)
    type(multifab) , intent(in   ) :: eta(:)              ! cell-centered
    type(multifab) , intent(in   ) :: eta_ed(:,:)         ! nodal (2D), edge-based (3D)
    type(multifab) , intent(in   ) :: temperature(:)      ! cell-centered
    type(multifab) , intent(in   ) :: temperature_ed(:,:) ! nodal (2D), edge-based (3D)
    real(dp_t)     , intent(in   ) :: dx(:,:),dt,weights(:)

    ! local
    integer :: n,nlevs,dm,i,comp
    integer :: ng_c,ng_e,ng_y,ng_w,ng_n,ng_f,ng_t,ng_m

    real(dp_t) :: variance

    type(multifab) :: mflux_cc_temp(mla%nlevel)
    type(multifab) :: mflux_nd_temp(mla%nlevel)
    type(multifab) :: mflux_ed_temp(mla%nlevel,3)

    logical :: nodal_temp(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:), dp(:,:,:,:), sp(:,:,:,:), tp(:,:,:,:)
    real(kind=dp_t), pointer :: fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:)
    real(kind=dp_t), pointer :: dxp(:,:,:,:), dyp(:,:,:,:), dzp(:,:,:,:)
    real(kind=dp_t), pointer :: ep1(:,:,:,:), ep2(:,:,:,:), ep3(:,:,:,:)
    real(kind=dp_t), pointer :: mp1(:,:,:,:), mp2(:,:,:,:), mp3(:,:,:,:)
    integer :: lo(mla%dim), hi(mla%dim)

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! create a place to store temporary copy of the random numbers
       ! these copies will be scaled and used to create the stochastic flux divergence
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! we need dm cell-centered fluxes for momentum
       call multifab_build(mflux_cc_temp(n),mla%la(n),dm,max(1,filtering_width))
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          nodal_temp = .true.
          call multifab_build(mflux_nd_temp(n),mla%la(n),2,filtering_width,nodal_temp)
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(mflux_ed_temp(n,1),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed_temp(n,2),mla%la(n),2,filtering_width,nodal_temp)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(mflux_ed_temp(n,3),mla%la(n),2,filtering_width,nodal_temp)
       end if

       call multifab_setval(mflux_cc_temp(n),0.d0,all=.true.)
       if (dm .eq. 2) then
          call multifab_setval(mflux_nd_temp(n),0.d0,all=.true.)
       else if (dm .eq. 3) then
          call multifab_setval(mflux_ed_temp(n,1),0.d0,all=.true.)
          call multifab_setval(mflux_ed_temp(n,2),0.d0,all=.true.)
          call multifab_setval(mflux_ed_temp(n,3),0.d0,all=.true.)
       end if

       ! add weighted contribution of fluxes
       do comp=1,n_rngs
          call saxpy(mflux_cc_temp(n),weights(comp),mflux_cc(n,comp),all=.true.)
          if (dm .eq. 2) then
             call saxpy(mflux_nd_temp(n),weights(comp),mflux_nd(n,comp),all=.true.)
          else if (dm .eq. 3) then
             call saxpy(mflux_ed_temp(n,1),weights(comp),mflux_ed(n,1,comp),all=.true.)
             call saxpy(mflux_ed_temp(n,2),weights(comp),mflux_ed(n,2,comp),all=.true.)
             call saxpy(mflux_ed_temp(n,3),weights(comp),mflux_ed(n,3,comp),all=.true.)
          end if
       end do

    end do

    ng_c = mflux_cc_temp(1)%ng
    ng_y = eta(1)%ng
    ng_t = temperature(1)%ng
    ng_f = stoch_m_force(1,1)%ng

    do n=1,nlevs

       ! include eta and temperature contribution in an ijk loop
       variance = sqrt(variance_coef_mom*2.d0*k_B/(product(dx(n,1:dm))*dt))

       if (dm .eq. 2) then

          ng_n = mflux_nd_temp(1)%ng
          ng_w = eta_ed(1,1)%ng
          ng_m = temperature_ed(1,1)%ng
          
          do i=1,nfabs(mflux_cc_temp(n))
             
             fp  => dataptr(mflux_cc_temp(n),i)
             sp  => dataptr(mflux_nd_temp(n),i)
             dp  => dataptr(eta(n),i)
             ep1 => dataptr(eta_ed(n,1),i)
             tp  => dataptr(temperature(n),i)
             mp1 => dataptr(temperature_ed(n,1),i)
             lo = lwb(get_box(mflux_cc_temp(n),i))
             hi = upb(get_box(mflux_cc_temp(n),i))
             ! multiply by variance
             fp = variance*fp
             sp = variance*sp
             ! multiply by sqrt(temperature*eta)
             call mult_by_sqrt_eta_2d(fp(:,:,1,:),ng_c,sp(:,:,1,:),ng_n, &
                                      dp(:,:,1,1),ng_y,ep1(:,:,1,1),ng_w, &
                                      tp(:,:,1,1),ng_t,mp1(:,:,1,1),ng_m,lo,hi)
             ! apply boundary conditions
             call mflux_bc_2d(sp(:,:,1,:),ng_n,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end do

          ! sync up random numbers at boundaries and ghost cells
          call multifab_internal_sync(mflux_nd_temp(n))
          call multifab_fill_boundary(mflux_nd_temp(n))
          call multifab_fill_boundary(mflux_cc_temp(n))
          
          if(filtering_width>0) then
             call multifab_filter(mflux_nd_temp(n), filtering_width)
             call multifab_filter(mflux_cc_temp(n), filtering_width)
             call multifab_fill_boundary(mflux_cc_temp(n)) ! First ghost cell is used in divergence
          end if

       else if (dm .eq. 3) then

          ng_e = mflux_ed_temp(1,1)%ng
          ng_w = eta_ed(1,1)%ng
          ng_m = temperature_ed(1,1)%ng

          do i=1,nfabs(mflux_cc_temp(n))
             fp  => dataptr(mflux_cc_temp(n),i)
             fxp => dataptr(mflux_ed_temp(n,1),i)
             fyp => dataptr(mflux_ed_temp(n,2),i)
             fzp => dataptr(mflux_ed_temp(n,3),i)
             dp => dataptr(eta(n),i)
             ep1 => dataptr(eta_ed(n,1),i)
             ep2 => dataptr(eta_ed(n,2),i)
             ep3 => dataptr(eta_ed(n,3),i)
             tp => dataptr(temperature(n),i)
             mp1 => dataptr(temperature_ed(n,1),i)
             mp2 => dataptr(temperature_ed(n,2),i)
             mp3 => dataptr(temperature_ed(n,3),i)
             lo = lwb(get_box(mflux_cc_temp(n),i))
             hi = upb(get_box(mflux_cc_temp(n),i))
             ! multiply by variance
             fp  = variance*fp
             fxp = variance*fxp
             fyp = variance*fyp
             fzp = variance*fzp
             ! multiply by sqrt(temperature*eta)
             call mult_by_sqrt_eta_3d(fp(:,:,:,:),ng_c, &
                                      fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e, &
                                      dp(:,:,:,1),ng_y, &
                                      ep1(:,:,:,1),ep2(:,:,:,1),ep3(:,:,:,1),ng_w, &
                                      tp(:,:,:,1),ng_t, &
                                      mp1(:,:,:,1),mp2(:,:,:,1),mp3(:,:,:,1),ng_m, &
                                      lo,hi)
             ! apply boundary conditions
             call mflux_bc_3d(fxp(:,:,:,:),fyp(:,:,:,:),fzp(:,:,:,:),ng_e,lo,hi, &
                              the_bc_level(n)%phys_bc_level_array(i,:,:))
          end do
          
          ! sync up random numbers at boundaries and ghost cells
          call multifab_internal_sync(mflux_ed_temp(n,1))
          call multifab_internal_sync(mflux_ed_temp(n,2))
          call multifab_internal_sync(mflux_ed_temp(n,3))
          call multifab_fill_boundary(mflux_ed_temp(n,1))
          call multifab_fill_boundary(mflux_ed_temp(n,2))
          call multifab_fill_boundary(mflux_ed_temp(n,3))
          call multifab_fill_boundary(mflux_cc_temp(n))

          if(filtering_width>0) then
             call multifab_filter(mflux_ed_temp(n,1), filtering_width)
             call multifab_filter(mflux_ed_temp(n,2), filtering_width)
             call multifab_filter(mflux_ed_temp(n,3), filtering_width)
             call multifab_filter(mflux_cc_temp(n)  , filtering_width)
             call multifab_fill_boundary(mflux_cc_temp(n)) ! First ghost cell is used in divergence
          end if

       end if

       ! calculate divergence and add to stoch_m_force
       do i=1,nfabs(stoch_m_force(n,1))
          fp => dataptr(mflux_cc_temp(n), i)
          dxp => dataptr(stoch_m_force(n,1),i)
          dyp => dataptr(stoch_m_force(n,2),i)
          lo =  lwb(get_box(stoch_m_force(n,1), i))
          hi =  upb(get_box(stoch_m_force(n,1), i))
          select case (dm)
          case (2)
             sp => dataptr(mflux_nd_temp(n), i)
             call stoch_m_force_2d(fp(:,:,1,:), sp(:,:,1,:), dxp(:,:,1,1), dyp(:,:,1,1), &
                                   ng_c, ng_n, ng_f, dx(n,:), lo, hi)
          case (3)
             dzp => dataptr(stoch_m_force(n,3), i)
             fxp => dataptr(mflux_ed_temp(n,1), i)
             fyp => dataptr(mflux_ed_temp(n,2), i)
             fzp => dataptr(mflux_ed_temp(n,3), i)
             call stoch_m_force_3d(fp(:,:,:,:), fxp(:,:,:,:), fyp(:,:,:,:), fzp(:,:,:,:), &
                                   dxp(:,:,:,1), dyp(:,:,:,1), dzp(:,:,:,1), &
                                   ng_c, ng_e, ng_f, dx(n,:), lo, hi)

          end select
       end do

       do i=1,dm
          call multifab_physbc_domainvel(stoch_m_force(n,i),vel_bc_comp+i-1, &
                                         the_bc_level(n),dx(n,:))
       end do

    end do

    do n=1,nlevs
       call multifab_destroy(mflux_cc_temp(n))
       if (dm .eq. 2) then
          call multifab_destroy(mflux_nd_temp(n))
       else if (dm .eq. 3) then
          call multifab_destroy(mflux_ed_temp(n,1))
          call multifab_destroy(mflux_ed_temp(n,2))
          call multifab_destroy(mflux_ed_temp(n,3))
       end if
    end do

  contains
    
    subroutine mult_by_sqrt_eta_2d(mflux_cc,ng_c,mflux_nd,ng_n,eta,ng_y,eta_nodal,ng_w, &
         temperature,ng_t,temperature_nodal,ng_m,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n,ng_y,ng_w,ng_t,ng_m
      real(kind=dp_t), intent(inout) ::          mflux_cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
      real(kind=dp_t), intent(inout) ::          mflux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      real(kind=dp_t), intent(in   ) ::               eta(lo(1)-ng_y:,lo(2)-ng_y:)
      real(kind=dp_t), intent(in   ) ::         eta_nodal(lo(1)-ng_w:,lo(2)-ng_w:)
      real(kind=dp_t), intent(in   ) ::       temperature(lo(1)-ng_t:,lo(2)-ng_t:)
      real(kind=dp_t), intent(in   ) :: temperature_nodal(lo(1)-ng_m:,lo(2)-ng_m:)

      ! local
      integer i,j

      do j=lo(2)-1,hi(2)+1
         do i=lo(1)-1,hi(1)+1
            mflux_cc(i,j,:) = mflux_cc(i,j,:) * sqrt(eta(i,j)*temperature(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            mflux_nd(i,j,:) = mflux_nd(i,j,:) * sqrt(eta_nodal(i,j)*temperature_nodal(i,j))
         end do
      end do

    end subroutine mult_by_sqrt_eta_2d

    subroutine mult_by_sqrt_eta_3d(mflux_cc,ng_c,mflux_xy,mflux_xz,mflux_yz,ng_e,eta,ng_y, &
                                   eta_xy,eta_xz,eta_yz,ng_w,temperature,ng_t, &
                                   temperature_xy,temperature_xz,temperature_yz,ng_m,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_e,ng_y,ng_w,ng_t,ng_m
      real(kind=dp_t), intent(inout) ::       mflux_cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
      real(kind=dp_t), intent(inout) ::       mflux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) ::       mflux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) ::       mflux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) ::            eta(lo(1)-ng_y:,lo(2)-ng_y:,lo(3)-ng_y:)
      real(kind=dp_t), intent(in   ) ::         eta_xy(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
      real(kind=dp_t), intent(in   ) ::         eta_xz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
      real(kind=dp_t), intent(in   ) ::         eta_yz(lo(1)-ng_w:,lo(2)-ng_w:,lo(3)-ng_w:)
      real(kind=dp_t), intent(in   ) ::    temperature(lo(1)-ng_t:,lo(2)-ng_t:,lo(3)-ng_t:)
      real(kind=dp_t), intent(in   ) :: temperature_xy(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(in   ) :: temperature_xz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(in   ) :: temperature_yz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)

      ! local
      integer i,j,k

      do k=lo(3)-1,hi(3)+1
         do j=lo(2)-1,hi(2)+1
            do i=lo(1)-1,hi(1)+1
               mflux_cc(i,j,k,:) = mflux_cc(i,j,k,:) * sqrt(eta(i,j,k)*temperature(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               mflux_xy(i,j,k,:) = mflux_xy(i,j,k,:) * sqrt(eta_xy(i,j,k)*temperature_xy(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               mflux_xz(i,j,k,:) = mflux_xz(i,j,k,:) * sqrt(eta_xz(i,j,k)*temperature_xz(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               mflux_yz(i,j,k,:) = mflux_yz(i,j,k,:) * sqrt(eta_yz(i,j,k)*temperature_yz(i,j,k))
            end do
         end do
      end do

    end subroutine mult_by_sqrt_eta_3d

    subroutine mflux_bc_2d(mflux_nd,ng_n,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_n
      real(kind=dp_t), intent(inout) :: mflux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      integer        , intent(in   ) :: phys_bc(:,:)

      ! y-mom fluxes that live on x-domain boundaries
      if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1),lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(lo(1),lo(2):hi(2)+1,:)
      else if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_RESERVOIR) then
         mflux_nd(lo(1),lo(2):hi(2)+1,:) = 0.d0
      end if
      if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = sqrt(2.d0)*mflux_nd(hi(1)+1,lo(2):hi(2)+1,:)
      else if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_RESERVOIR) then
         mflux_nd(hi(1)+1,lo(2):hi(2)+1,:) = 0.d0
      end if

      ! x-mom fluxes that live on y-domain boundaries
      if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,lo(2),:) = sqrt(2.d0)*mflux_nd(lo(1):hi(1)+1,lo(2),:)
      else if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,lo(2),:) = 0.d0
      end if
      if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,hi(2)+1,:) = sqrt(2.d0)*mflux_nd(lo(1):hi(1)+1,hi(2)+1,:)
      else if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_RESERVOIR) then
         mflux_nd(lo(1):hi(1)+1,hi(2)+1,:) = 0.d0
      end if

    end subroutine mflux_bc_2d

    subroutine mflux_bc_3d(mflux_xy,mflux_xz,mflux_yz,ng_e,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_e
      real(kind=dp_t), intent(inout) :: mflux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) :: mflux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      integer        , intent(in   ) :: phys_bc(:,:)

      ! y-mom and z-mom fluxes that live on x-domain boundaries
      if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:)
         mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:)
      else if (phys_bc(1,1) .eq. SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_RESERVOIR) then
         mflux_xy(lo(1),lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_xz(lo(1),lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
      end if
      if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:)
         mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:)
      else if (phys_bc(1,2) .eq. SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_RESERVOIR) then
         mflux_xy(hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_xz(hi(1)+1,lo(2):hi(2),lo(3):hi(3)+1,:) = 0.d0
      end if

      ! x-mom and z-mom fluxes that live on y-domain boundaries
      if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:)
         mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:)
      else if (phys_bc(2,1) .eq. SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,lo(2),lo(3):hi(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2),lo(3):hi(3)+1,:) = 0.d0
      end if
      if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = sqrt(2.d0)*mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:)
         mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:)
      else if (phys_bc(2,2) .eq. SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_RESERVOIR) then
         mflux_xy(lo(1):hi(1)+1,hi(2)+1,lo(3):hi(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),hi(2)+1,lo(3):hi(3)+1,:) = 0.d0
      end if

      ! x-mom and y-mom fluxes that live on z-domain boundaries
      if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = sqrt(2.d0)*mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:)
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:)
      else if (phys_bc(3,1) .eq. SLIP_WALL .or. phys_bc(3,1) .eq. SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),lo(3),:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,lo(3),:) = 0.d0
      end if
      if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:) = sqrt(2.d0)*mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:)
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:) = sqrt(2.d0)*mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:)
      else if (phys_bc(3,2) .eq. SLIP_WALL .or. phys_bc(3,2) .eq. SLIP_RESERVOIR) then
         mflux_xz(lo(1):hi(1)+1,lo(2):hi(2),hi(3)+1,:) = 0.d0
         mflux_yz(lo(1):hi(1),lo(2):hi(2)+1,hi(3)+1,:) = 0.d0
      end if

    end subroutine mflux_bc_3d

    subroutine stoch_m_force_2d(flux_cc,flux_nd,divx,divy,ng_c,ng_n,ng_f,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_n,ng_f
      real(kind=dp_t), intent(in   ) :: flux_cc(lo(1)-ng_c:,lo(2)-ng_c:,:)
      real(kind=dp_t), intent(in   ) :: flux_nd(lo(1)-ng_n:,lo(2)-ng_n:,:)
      real(kind=dp_t), intent(inout) ::    divx(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divy(lo(1)-ng_f:,lo(2)-ng_f:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      integer :: i,j
      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      ! divergence on x-faces
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1

            divx(i,j) = divx(i,j) + &
                 (flux_cc(i,j,1) - flux_cc(i-1,j,1)) * dxinv + &
                 (flux_nd(i,j+1,1) - flux_nd(i,j,1)) * dxinv

         end do
      end do

      ! divergence on y-faces
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

            divy(i,j) = divy(i,j) + &
                 (flux_nd(i+1,j,2) - flux_nd(i,j,2)) * dxinv + &
                 (flux_cc(i,j,2) - flux_cc(i,j-1,2)) * dxinv

         end do
      end do


    end subroutine stoch_m_force_2d

    subroutine stoch_m_force_3d(flux_cc,flux_xy,flux_xz,flux_yz,divx,divy,divz, &
                                ng_c,ng_e,ng_f,dx,lo,hi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_c,ng_e,ng_f
      real(kind=dp_t), intent(in   ) :: flux_cc(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:,:)
      real(kind=dp_t), intent(in   ) :: flux_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) :: flux_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(in   ) :: flux_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
      real(kind=dp_t), intent(inout) ::    divx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(inout) ::    divz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k
      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      ! divergence on x-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1

               divx(i,j,k) = divx(i,j,k) + &
                    (flux_cc(i,j,k,1) - flux_cc(i-1,j,k,1)) * dxinv + &
                    (flux_xy(i,j+1,k,1) - flux_xy(i,j,k,1)) * dxinv + &
                    (flux_xz(i,j,k+1,1) - flux_xz(i,j,k,1)) * dxinv

            end do
         end do
      end do

      ! divergence on y-faces
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)

               divy(i,j,k) = divy(i,j,k) + &
                    (flux_xy(i+1,j,k,2) - flux_xy(i,j,k,2)) * dxinv + &
                    (flux_cc(i,j,k,2) - flux_cc(i,j-1,k,2)) * dxinv + &
                    (flux_yz(i,j,k+1,1) - flux_yz(i,j,k,1)) * dxinv

            end do
         end do
      end do

      ! divergence on z-faces
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)

               divz(i,j,k) = divz(i,j,k) + &
                    (flux_xz(i+1,j,k,2) - flux_xz(i,j,k,2)) * dxinv + &
                    (flux_yz(i,j+1,k,2) - flux_yz(i,j,k,2)) * dxinv + &
                    (flux_cc(i,j,k,3) - flux_cc(i,j,k-1,3)) * dxinv

            end do
         end do
      end do

    end subroutine stoch_m_force_3d

  end subroutine stochastic_m_fluxdiv

  ! fill the stochastic multifabs with random numbers
  subroutine fill_m_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,dm,comp

    nlevs = mla%nlevel
    dm = mla%dim

    do comp=1,n_rngs

       ! Diagonal components of stochastic stress tensor for momentum
       select case(stoch_stress_form)
       case (0) ! Non-symmetric
          call multifab_fill_random(mflux_cc(:,comp))
       case default ! Symmetric
          call multifab_fill_random(mflux_cc(:,comp), variance=2.d0)
       end select

       ! Off-diagonal components of stochastic stress tensor for momentum
       if (dm .eq. 2) then
          ! in 2D, we need 2 random fluxes at each node
          select case(stoch_stress_form)
          case(0) ! Non-symmetric
             call multifab_fill_random(mflux_nd(:,comp))
          case default ! Symmetric
             call multifab_fill_random(mflux_nd(:,comp), comp=1)
             do n=1,nlevs
                call multifab_copy_c(mflux_nd(n,comp),2,mflux_nd(n,comp),1)
             end do
          end select
       else if (dm .eq. 3) then
          ! in 3D, we need 2 random fluxes at each edge
          select case(stoch_stress_form)
          case(0) ! Non-symmetric
             call multifab_fill_random(mflux_ed(:,1,comp))
             call multifab_fill_random(mflux_ed(:,2,comp))
             call multifab_fill_random(mflux_ed(:,3,comp))
          case default ! Symmetric
             call multifab_fill_random(mflux_ed(:,1,comp), comp=1)
             call multifab_fill_random(mflux_ed(:,2,comp), comp=1)
             call multifab_fill_random(mflux_ed(:,3,comp), comp=1)
             do n = 1, nlevs
                call multifab_copy_c(mflux_ed(n,1,comp),2,mflux_ed(n,1,comp),1)
                call multifab_copy_c(mflux_ed(n,2,comp),2,mflux_ed(n,2,comp),1)
                call multifab_copy_c(mflux_ed(n,3,comp),2,mflux_ed(n,3,comp),1)
             end do
          end select
       end if

    end do

  end subroutine fill_m_stochastic

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_m_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,dm,comp
    logical :: nodal_temp(mla%dim)
    
    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(mflux_cc(mla%nlevel         , n_rngs))
    allocate(mflux_nd(mla%nlevel         , n_rngs))
    allocate(mflux_ed(mla%nlevel, 3      , n_rngs))

    do n=1,nlevs
       do comp=1,n_rngs
          ! we need dm cell-centered fluxes for momentum
          call multifab_build(mflux_cc(n,comp),mla%la(n),dm,max(1,filtering_width))
          if (dm .eq. 2) then
             ! in 2D, we need 2 random fluxes at each node
             nodal_temp = .true.
             call multifab_build(mflux_nd(n,comp),mla%la(n),2,filtering_width,nodal_temp)
          else if (dm .eq. 3) then
             ! in 3D, we need 2 random fluxes at each edge
             nodal_temp(1) = .true.
             nodal_temp(2) = .true.
             nodal_temp(3) = .false.
             call multifab_build(mflux_ed(n,1,comp),mla%la(n),2,filtering_width,nodal_temp)
             nodal_temp(1) = .true.
             nodal_temp(2) = .false.
             nodal_temp(3) = .true.
             call multifab_build(mflux_ed(n,2,comp),mla%la(n),2,filtering_width,nodal_temp)
             nodal_temp(1) = .false.
             nodal_temp(2) = .true.
             nodal_temp(3) = .true.
             call multifab_build(mflux_ed(n,3,comp),mla%la(n),2,filtering_width,nodal_temp)
          end if
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

  end subroutine init_m_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_m_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,dm,comp
    
    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do comp=1,n_rngs
          call multifab_destroy(mflux_cc(n,comp))
          if (dm .eq. 2) then
             call multifab_destroy(mflux_nd(n,comp))
          else if (dm .eq. 3) then
             call multifab_destroy(mflux_ed(n,1,comp))
             call multifab_destroy(mflux_ed(n,2,comp))
             call multifab_destroy(mflux_ed(n,3,comp))
          end if
       end do
    end do
    
    deallocate(mflux_cc,mflux_nd,mflux_ed)

  end subroutine destroy_m_stochastic

 ! Add equilibrium fluctuations to the momentum (valid and ghost regions)
 subroutine add_m_fluctuations_1(mla,dx,variance,s_face,temperature_face,m_face)

   type(ml_layout), intent(in   ) :: mla
   real(dp_t)     , intent(in   ) :: variance, dx(:,:)
   type(multifab) , intent(in   ) :: s_face(:,:), temperature_face(:,:)
   type(multifab) , intent(inout) :: m_face(:,:)

   ! local
   type(multifab) :: mactemp(mla%nlevel,mla%dim)
   integer :: n,i,dm,nlevs
   real(dp_t) :: av_mom(mla%dim)

   dm = mla%dim
   nlevs = mla%nlevel

   do n=1,nlevs
      do i=1,dm
         call multifab_build_edge(mactemp(n,i),mla%la(n),1,1,i)
      end do
   end do

   ! Generate random numbers first and store them in mactemp
   do n=1,nlevs
      do i=1,dm
         call multifab_fill_random(mactemp(n:n,i), &
                                   variance=abs(variance)*k_B/product(dx(n,1:dm)), &
                                   variance_mfab=s_face(n:n,i), &
                                   variance_mfab2=temperature_face(n:n,i))
         call saxpy(m_face(n,i), 1.0_dp_t, mactemp(n,i), all=.true.)
      end do
   end do

   do n=1,nlevs
      do i=1,dm
         ! We need to ensure periodic BCs are obeyed for the random values
         call multifab_internal_sync(m_face(n,i))
      end do
   enddo
   
   if(variance<0) then ! Ensure zero total momentum
      if (parallel_IOProcessor()) then
         write(*,"(A,100G17.9)") "Randomly INITIALized momenta"
      end if
      
      call sum_momenta(mla, m_face, av_mom)
      do i=1,dm
         call setval(mactemp(1,i), -av_mom(i))
         call saxpy(m_face(1,i), 1.0_dp_t, mactemp(1,i))
      end do
   end if

   do n=1,nlevs
      do i=1,dm
         call multifab_destroy(mactemp(n,i))
      end do
   end do

 end subroutine add_m_fluctuations_1

 subroutine add_m_fluctuations_2(mla,dx,variance,umac,rhotot,Temp,the_bc_tower)

   type(ml_layout), intent(in   ) :: mla
   real(dp_t)     , intent(in   ) :: variance, dx(:,:)
   type(multifab) , intent(inout) :: umac(:,:)
   type(multifab) , intent(in   ) :: rhotot(:)
   type(multifab) , intent(in   ) :: Temp(:)
   type(bc_tower) , intent(in   ) :: the_bc_tower

   ! local
   type(multifab) ::      mold(mla%nlevel,mla%dim)
   type(multifab) :: rhotot_fc(mla%nlevel,mla%dim)
   type(multifab) ::   Temp_fc(mla%nlevel,mla%dim)

   integer :: i,dm,n,nlevs

   nlevs = mla%nlevel
   dm = mla%dim
      
   ! temporary multifabs
   do n=1,nlevs
      do i=1,dm
         call multifab_build_edge(mold(n,i),mla%la(n),1,1,i)
         call multifab_build_edge(rhotot_fc(n,i),mla%la(n),1,1,i)
         call multifab_build_edge(Temp_fc(n,i),mla%la(n),1,1,i)
      end do
   end do

   ! compute rhotot on faces
   call average_cc_to_face(nlevs,rhotot,rhotot_fc,1,scal_bc_comp,1, &
                           the_bc_tower%bc_tower_array)

   ! compute Temp on faces
   call average_cc_to_face(nlevs,Temp,Temp_fc,1,temp_bc_comp,1, &
                           the_bc_tower%bc_tower_array)

   ! compute mold
   call convert_m_to_umac(mla,rhotot_fc,mold,umac,.false.)

   ! add fluctuations to mold and convert back to umac
   call add_m_fluctuations(mla,dx,variance,rhotot_fc,Temp_fc,mold)

   ! convert back to umac
   call convert_m_to_umac(mla,rhotot_fc,mold,umac,.true.)

   do n=1,nlevs
      do i=1,dm
         call multifab_destroy(mold(n,i))
         call multifab_destroy(rhotot_fc(n,i))
         call multifab_destroy(Temp_fc(n,i))
      end do
   end do

 end subroutine add_m_fluctuations_2

end module stochastic_m_fluxdiv_module
