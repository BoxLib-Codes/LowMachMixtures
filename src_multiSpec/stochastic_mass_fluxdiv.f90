module stochastic_mass_fluxdiv_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use mass_flux_utilities_module
  use multifab_fill_random_module
  use div_and_grad_module
  use convert_stag_module
  use matvec_mul_module
  use correction_flux_module
  use zero_edgeval_module
  use probin_common_module, only: k_B, variance_coef_mass
  use probin_multispecies_module, only: nspecies, correct_flux, is_nonisothermal

  implicit none

  private

  public :: stochastic_mass_fluxdiv, fill_mass_stochastic, init_mass_stochastic, destroy_mass_stochastic

  ! stochastic fluxes for mass densities are face-centered
  type(multifab), allocatable, save :: stoch_W_fc(:,:,:)

  integer, save :: n_rngs ! how many random number stages
  
contains
  
  subroutine stochastic_mass_fluxdiv(mla,rho,rhotot,molarconc,molmtot,chi,&
                                     Gama,stoch_fluxdiv,flux_total,dx,dt,weights, &
                                     the_bc_level,increment_in)

    type(ml_layout), intent(in   )   :: mla
    type(multifab) , intent(in   )   :: rho(:)
    type(multifab) , intent(in   )   :: rhotot(:)
    type(multifab) , intent(in   )   :: molarconc(:)
    type(multifab) , intent(in   )   :: molmtot(:)
    type(multifab) , intent(in   )   :: chi(:)
    type(multifab) , intent(in   )   :: Gama(:)
    type(multifab) , intent(inout)   :: stoch_fluxdiv(:)
    type(multifab) , intent(inout)   :: flux_total(:,:)
    real(kind=dp_t), intent(in   )   :: dx(:,:)
    real(kind=dp_t), intent(in   )   :: dt
    real(kind=dp_t), intent(in   )   :: weights(:)         
    type(bc_level) , intent(in   )   :: the_bc_level(:)
    logical  ,  intent(in), optional :: increment_in

    ! Local variables
    type(multifab)   :: Lonsager(mla%nlevel)            ! cholesky factored Lonsager 
    type(multifab)   :: Lonsager_fc(mla%nlevel,mla%dim) ! cholesky factored Lonsager on face
    type(multifab)   :: flux(mla%nlevel,mla%dim)        ! face-centered stochastic flux
    integer          :: n,nlevs,i,dm,rng
    real(kind=dp_t)  :: variance
    logical :: increment

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stochastic_mass_fluxdiv")

    ! do we increment or overwrite div?
    increment = .false.
    if (present(increment_in)) increment = increment_in

    nlevs = mla%nlevel
    dm    = mla%dim
 
    ! populate the variance (only first level)
    variance = sqrt(2.d0*k_B*variance_coef_mass/(product(dx(1,1:dm))*dt))

    ! build multifabs 
    do n=1,nlevs
       call multifab_build(Lonsager(n), mla%la(n), nspecies**2, rho(n)%ng)
       do i=1,dm
          call multifab_build_edge(Lonsager_fc(n,i),   mla%la(n), nspecies**2, 0, i)
          call multifab_build_edge(flux(n,i), mla%la(n), nspecies,    0, i)
       end do
    end do
 
    ! set flux to zero
    do n=1,nlevs
       do i = 1,dm
          call setval(flux(n,i), 0.d0, all=.true.)   
       end do   
    end do   
    
    ! convert stoch_W_fc into flux
    do n=1,nlevs
       do i = 1,dm
          do rng=1,size(weights)
             call saxpy(flux(n,i), weights(rng), stoch_W_fc(n,i,rng))
          end do   
       end do   
    end do
    
    ! compute cell-centered cholesky-factored Lonsager^(1/2)
    call compute_Lonsager(mla,rho,rhotot,molarconc,molmtot,chi,Gama,Lonsager)
                  
    ! compute face-centered cholesky factor of cell-centered cholesky factored Lonsager^(1/2)
    call average_cc_to_face(nlevs,Lonsager,Lonsager_fc,1,tran_bc_comp,nspecies**2,the_bc_level,.false.)

    ! compute variance X cholesky-Lonsager-face X W(0,1) 
    do n=1,nlevs
       do i=1,dm
          call matvec_mul(mla, flux(n,i), Lonsager_fc(n,i), nspecies)
          call multifab_mult_mult_s(flux(n,i), variance, 0)
       end do
    end do  
    
    ! sync the fluxes at the boundaries
    do n=1,nlevs
       do i=1,dm
          call multifab_internal_sync(flux(n,i))
          call multifab_fill_boundary(flux(n,i))  
       end do
    end do

    ! If there are walls with zero-flux boundary conditions
    if(is_nonisothermal) then
       do n=1,nlevs
          call zero_edgeval_walls(flux(n,:),1,nspecies,the_bc_level(n))
       end do   
    end if

    !correct fluxes to ensure mass conservation to roundoff
    if (correct_flux .and. (nspecies .gt. 1)) then
       !write(*,*) "Checking conservation of stochastic fluxes"
       call correction_flux(mla, rho, rhotot, flux, the_bc_level)
    end if

    ! add fluxes to flux_total
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(flux_total(n,i),1,flux(n,i),1,nspecies,0)
       end do
    end do
 
    ! compute divergence of stochastic flux
    call compute_div(mla,flux,stoch_fluxdiv,dx,1,1,nspecies,increment)

    ! free the multifab allocated memory
    do n=1,nlevs
       call multifab_destroy(Lonsager(n))
       do i=1,dm
          call multifab_destroy(Lonsager_fc(n,i))
          call multifab_destroy(flux(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine stochastic_mass_fluxdiv

  ! call this once at the beginning of simulation to allocate multifabs
  ! that will hold random numbers
  subroutine init_mass_stochastic(mla,n_rngs_in)

    type(ml_layout), intent(in   ) :: mla
    integer        , intent(in   ) :: n_rngs_in

    ! local
    integer :: n,nlevs,i,dm,comp
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"init_mass_stochastic")

    n_rngs = n_rngs_in

    nlevs = mla%nlevel
    dm = mla%dim

    allocate(stoch_W_fc(mla%nlevel, mla%dim, n_rngs))

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             ! we need one face-centered flux for each concentration
             call multifab_build_edge(stoch_W_fc(n,i,comp),mla%la(n),nspecies,0,i)
          end do
       end do ! end loop over n_rngs
    end do ! end loop over nlevs

    call destroy(bpt)

  end subroutine init_mass_stochastic

  ! call this once at the end of simulation to deallocate memory
  subroutine destroy_mass_stochastic(mla)

    type(ml_layout), intent(in   ) :: mla

    ! local
    integer :: n,nlevs,i,dm,comp
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"destroy_mass_stochastic")

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do comp=1,n_rngs
          do i=1,dm
             call multifab_destroy(stoch_W_fc(n,i,comp))
          end do
       end do
    end do
    
    deallocate(stoch_W_fc)

    call destroy(bpt)

  end subroutine destroy_mass_stochastic

  subroutine fill_mass_stochastic(mla,the_bc_level)
  
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! Local variables
    integer :: dm,nlevs,box,i,rng
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_mass_stochastic")

    nlevs = mla%nlevel
    dm    = mla%dim    
    
    ! generate and store the stochastic flux (random numbers)
    do rng=1, n_rngs
       do i = 1,dm
          call multifab_fill_random(stoch_W_fc(:,i,rng))
       end do
    end do

    ! apply boundary conditions to stochastic fluxes
    call stoch_mass_bc(mla,the_bc_level)
  
    call destroy(bpt)

  end subroutine fill_mass_stochastic

  subroutine stoch_mass_bc(mla,the_bc_level)
    
    type(ml_layout), intent(in   )  :: mla
    type(bc_level) , intent(in   )  :: the_bc_level(:)

    ! local
    integer :: n,nlevs,dm,idim,comp,i,ng_f
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"stoch_mass_bc")

    nlevs = mla%nlevel
    dm = mla%dim
    
    ng_f = stoch_W_fc(1,1,1)%ng

    do n=1,nlevs
       do idim=1,dm
          do comp=1,n_rngs
             do i=1,nfabs(stoch_W_fc(n,idim,comp))
                fp => dataptr(stoch_W_fc(n,idim,comp),i)
                lo = lwb(get_box(stoch_W_fc(n,idim,comp),i))
                hi = upb(get_box(stoch_W_fc(n,idim,comp),i))
                select case (dm)
                case (2)
                   call stoch_mass_bc_2d(fp(:,:,1,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                case (3)
                   call stoch_mass_bc_3d(fp(:,:,:,:),ng_f,idim,lo,hi, &
                                         the_bc_level(n)%phys_bc_level_array(i,:,:))
                end select
             end do
          end do
       end do
    end do

    call destroy(bpt)

  end subroutine stoch_mass_bc

  subroutine stoch_mass_bc_2d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),:)
         end if

      else

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_2d

  subroutine stoch_mass_bc_3d(sflux,ng_f,idim,lo,hi,phys_bc)

      integer        , intent(in   ) :: lo(:),hi(:),ng_f,idim
      integer        , intent(in   ) :: phys_bc(:,:)
      real(kind=dp_t), intent(inout) :: sflux(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:,:)

      if (idim .eq. 1) then

         if (phys_bc(1,1) .eq. NO_SLIP_WALL .or. phys_bc(1,1) .eq. SLIP_WALL) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_WALL .or. phys_bc(1,2) .eq. SLIP_WALL) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(1,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1),lo(2):hi(2),lo(3):hi(3),:)
         end if

         if (phys_bc(1,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(hi(1)+1,lo(2):hi(2),lo(3):hi(3),:)
         end if

      else if (idim .eq. 2) then

         if (phys_bc(2,1) .eq. NO_SLIP_WALL .or. phys_bc(2,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_WALL .or. phys_bc(2,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = 0.d0
         end if

         if (phys_bc(2,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2),lo(3):hi(3),:)
         end if

         if (phys_bc(2,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),hi(2)+1,lo(3):hi(3),:)
         end if

      else

         if (phys_bc(3,1) .eq. NO_SLIP_WALL .or. phys_bc(3,1) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = 0.d0
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_WALL .or. phys_bc(3,2) .eq. SLIP_WALL) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = 0.d0
         end if

         if (phys_bc(3,1) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),lo(3),:)
         end if

         if (phys_bc(3,2) .eq. NO_SLIP_RESERVOIR) then
            sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:) = sqrt(2.0d0)*sflux(lo(1):hi(1),lo(2):hi(2),hi(3)+1,:)
         end if

      end if

  end subroutine stoch_mass_bc_3d
   
end module stochastic_mass_fluxdiv_module
