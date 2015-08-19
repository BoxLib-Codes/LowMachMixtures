module advance_timestep_inertial_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use convert_stag_module
  use mk_advective_s_fluxdiv_module
  use diffusive_m_fluxdiv_module
  use stochastic_m_fluxdiv_module
  use stochastic_mass_fluxdiv_module
  use compute_mass_fluxdiv_module
  use compute_HSE_pres_module
  use convert_m_to_umac_module
  use convert_rhoc_to_c_module
  use mk_advective_m_fluxdiv_module
  use reservoir_bc_fill_module
  use bds_module
  use gmres_module
  use div_and_grad_module
  use eos_check_module
  use mk_grav_force_module
  use compute_mixture_properties_module
  use mass_flux_utilities_module
  use multifab_physbc_module
  use multifab_physbc_extrap_module
  use multifab_physbc_stag_module
  use zero_edgeval_module
  use fill_rho_ghost_cells_module
  use probin_common_module, only: advection_type, grav, rhobar, variance_coef_mass, &
                                  variance_coef_mom, restart, barodiffusion_type
  use probin_gmres_module, only: gmres_abs_tol, gmres_rel_tol
  use probin_multispecies_module, only: nspecies
  use analysis_module

  implicit none

  private

  public :: advance_timestep_inertial

  ! special inhomogeneous boundary condition multifab
  ! vel_bc_n(nlevs,dm) are the normal velocities
  ! in 2D, vel_bc_t(nlevs,2) respresents
  !   1. y-velocity bc on x-faces (nodal)
  !   2. x-velocity bc on y-faces (nodal)
  ! in 3D, vel_bc_t(nlevs,6) represents
  !   1. y-velocity bc on x-faces (nodal in y and x)
  !   2. z-velocity bc on x-faces (nodal in z and x)
  !   3. x-velocity bc on y-faces (nodal in x and y)
  !   4. z-velocity bc on y-faces (nodal in z and y)
  !   5. x-velocity bc on z-faces (nodal in x and z)
  !   6. y-velocity bc on z-faces (nodal in y and z)
  type(multifab), allocatable, save :: vel_bc_n(:,:)
  type(multifab), allocatable, save :: vel_bc_t(:,:)

contains

  subroutine advance_timestep_inertial(mla,umac,rho_old,rho_new,rhotot_old,rhotot_new, &
                                       gradp_baro,pi,eta,eta_ed,kappa,Temp,Temp_ed, &
                                       diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                       dx,dt,time,the_bc_tower,istep)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho_old(:)
    type(multifab) , intent(inout) :: rho_new(:)
    type(multifab) , intent(inout) :: rhotot_old(:)
    type(multifab) , intent(inout) :: rhotot_new(:)
    type(multifab) , intent(inout) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: pi(:)
    ! eta and kappa need to enter consistent with old and leave consistent with new
    type(multifab) , intent(inout) :: eta(:)
    type(multifab) , intent(inout) :: eta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: kappa(:)
    type(multifab) , intent(inout) :: Temp(:)
    type(multifab) , intent(inout) :: Temp_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt,time
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ) :: istep

    ! local
    type(multifab) ::  rho_update(mla%nlevel)
    type(multifab) ::   bds_force(mla%nlevel)
    type(multifab) :: gmres_rhs_p(mla%nlevel)
    type(multifab) ::         dpi(mla%nlevel)
    type(multifab) ::        divu(mla%nlevel)
    type(multifab) ::        conc(mla%nlevel)
    type(multifab) ::  rho_nd_old(mla%nlevel)
    type(multifab) ::     rho_tmp(mla%nlevel)
    type(multifab) ::      p_baro(mla%nlevel)

    type(multifab) ::          mold(mla%nlevel,mla%dim)
    type(multifab) ::         mtemp(mla%nlevel,mla%dim)
    type(multifab) ::   m_a_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::   m_d_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::   m_s_fluxdiv(mla%nlevel,mla%dim)
    type(multifab) ::   gmres_rhs_v(mla%nlevel,mla%dim)
    type(multifab) ::         dumac(mla%nlevel,mla%dim)
    type(multifab) ::      umac_tmp(mla%nlevel,mla%dim)
    type(multifab) :: rhotot_fc_old(mla%nlevel,mla%dim)
    type(multifab) ::        gradpi(mla%nlevel,mla%dim)
    type(multifab) ::        rho_fc(mla%nlevel,mla%dim)
    type(multifab) ::     rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) ::    flux_total(mla%nlevel,mla%dim)

    integer :: i,dm,n,nlevs

    real(kind=dp_t) :: theta_alpha, norm_pre_rhs, gmres_abs_tol_in

    real(kind=dp_t) :: weights(1)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep_inertial")

    weights(1) = 1.d0

    nlevs = mla%nlevel
    dm = mla%dim

    theta_alpha = 1.d0/dt
    
    call build_bc_multifabs(mla)
    
    do n=1,nlevs
       call multifab_build( rho_update(n),mla%la(n),nspecies,0)
       call multifab_build(  bds_force(n),mla%la(n),nspecies,1)
       call multifab_build(gmres_rhs_p(n),mla%la(n),1       ,0)
       call multifab_build(         dpi(n),mla%la(n),1       ,1)
       call multifab_build(       divu(n),mla%la(n),1       ,0)
       call multifab_build(       conc(n),mla%la(n),nspecies,rho_old(n)%ng)
       call multifab_build(     p_baro(n),mla%la(n),1       ,1)
       do i=1,dm
          call multifab_build_edge(         mold(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(        mtemp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(  m_a_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  m_d_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  m_s_fluxdiv(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  gmres_rhs_v(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(        dumac(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(     umac_tmp(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(       gradpi(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(rhotot_fc_old(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(       rho_fc(n,i),mla%la(n),nspecies,0,i)
          call multifab_build_edge(    rhotot_fc(n,i),mla%la(n),1       ,1,i)
          call multifab_build_edge(   flux_total(n,i),mla%la(n),nspecies,0,i)
       end do
    end do

    ! make copies of old quantities
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(umac_tmp(n,i),1,umac(n,i),1,1,1)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 1 - Calculate Predictor Diffusive and Stochastic Fluxes
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! diff/stoch_mass_fluxdiv already contain F_i
    ! this was already done in Step 0 (initialization) or Step 6 from the previous time step

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 2 - Predictor Euler Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! average rho_old and rhotot_old to faces
    call average_cc_to_face(nlevs,   rho_old,   rho_fc    ,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_old,rhotot_fc_old,1,    scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! add D^n and St^n to rho_update
    do n=1,nlevs
       call setval(rho_update(n),0.d0,all=.true.)
       call multifab_plus_plus_c(rho_update(n),1, diff_mass_fluxdiv(n),1,nspecies,0)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies,0)
       end if
    end do

    ! add A^n to rho_update
    if (advection_type .ge. 1) then
      do n=1,nlevs
         ! set to zero to make sure ghost cells behind physical boundaries don't have NaNs
         call setval(bds_force(n),0.d0,all=.true.)
         call multifab_copy_c(bds_force(n),1,rho_update(n),1,nspecies,0)
         call multifab_fill_boundary(bds_force(n))
      end do

      if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          ! rho_fc (computed above) and rho_nd_old (computed here) are used to set boundary conditions
          do n=1,nlevs
             call multifab_build_nodal(rho_nd_old(n),mla%la(n),nspecies,1)
          end do
          call average_cc_to_node(nlevs,rho_old,rho_nd_old,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)

          ! the input s_tmp needs to have ghost cells filled with multifab_physbc_extrap
          ! instead of multifab_physbc
          do n=1,nlevs
             call multifab_build(rho_tmp(n),mla%la(n),nspecies,rho_old(n)%ng)
             call multifab_copy(rho_tmp(n),rho_old(n),rho_tmp(n)%ng)
             call multifab_physbc_extrap(rho_tmp(n),1,c_bc_comp,nspecies, &
                                         the_bc_tower%bc_tower_array(n))
          end do

          call bds(mla,umac,rho_tmp,rho_update,bds_force,rho_fc,rho_nd_old,dx,dt,1,nspecies, &
                   c_bc_comp,the_bc_tower,proj_type_in=2)

      else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=2)
      end if
    else
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,dx,1,nspecies)
    end if

    ! set rho_new = rho_old + dt * (A^n + D^n + St^n)
    do n=1,nlevs
       call multifab_mult_mult_s_c(rho_update(n),1,dt,nspecies,0)
       call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
       call multifab_plus_plus_c(rho_new(n),1,rho_update(n),1,nspecies,0)
    end do

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! rho to conc - NO GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)
    
    do n=1,nlevs
       call fill_rho_ghost_cells(conc(n),rhotot_new(n),the_bc_tower%bc_tower_array(n))
    end do

    ! conc to rho - INCLUDING GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.false.)

    ! average rho_new and rhotot_new to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,    scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 3 - Calculate Corrector Diffusive and Stochastic Fluxes
    ! Step 4 - Predictor Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt

    ! compute mold
    call convert_m_to_umac(mla,rhotot_fc_old,mold,umac,.false.)

    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute grad pi^n
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    if (barodiffusion_type .eq. 2) then
       ! barodiffusion uses lagged grad(pi)
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do
    else if (barodiffusion_type .eq. 3) then
       ! compute p0 from rho0*g
       call compute_HSE_pres(mla,rhotot_new,p_baro,dx,the_bc_tower)
       call compute_grad(mla,p_baro,gradp_baro,dx,1,pres_bc_comp,1,1, &
                         the_bc_tower%bc_tower_array)
    end if

    ! subtract grad pi^n from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! compute m_a_fluxdiv = A^n for momentum
    do n=1,nlevs
       do i=1,dm
          call setval(m_a_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do
    call mk_advective_m_fluxdiv(mla,umac,mold,m_a_fluxdiv,dx, &
                                the_bc_tower%bc_tower_array)

    ! add A^n for momentum to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_d_fluxdiv = A_0^n v^n
    do n=1,nlevs
       do i=1,dm
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute m_s_fluxdiv = div(Sigma^n)
    do n=1,nlevs
       do i=1,dm
          call setval(m_s_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,eta_ed, &
                                 Temp,Temp_ed,dx,dt,weights)
    end if

    ! add div(Sigma^n) to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add rho^n*g to gmres_rhs_v
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,rhotot_fc_old,rhotot_fc_old,the_bc_tower)
    end if

    ! compute (eta,kappa)^{*,n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time+dt, &
                                   the_bc_tower%bc_tower_array)

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F" so we later multiply by -1
    call compute_mass_fluxdiv(mla,rho_new,gradp_baro, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                      Temp,flux_total,dt,time,dx,weights, &
                                      the_bc_tower)

    do n=1,nlevs
       call multifab_mult_mult_s_c(diff_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       end if
       do i=1,dm
          call multifab_mult_mult_s_c(flux_total(n,i),1,-1.d0,nspecies,0)
       end do
    end do

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_tower%bc_tower_array)

    ! compute gmres_rhs_p
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
       end do
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^n vbar^n to gmres_rhs_v and add div vbar^n to gmres_rhs_p we
    ! are automatically putting the system in delta form WITH homogeneous boundary conditions
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

   ! compute mtemp = rho^{*,n+1} * vbar^n
   call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

   do n=1,nlevs
      do i=1,dm

         ! multiply mtemp by 1/dt
         call multifab_mult_mult_s_c(mtemp(n,i),1,1.d0/dt,1,0)

         ! subtract rho^{*,n+1} * vbar^n / dt from gmres_rhs_v
         call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)

      end do
   end do

   do n=1,nlevs
      do i=1,dm
         ! reset mtemp
         call multifab_setval(mtemp(n,i),0.d0,all=.true.)
      end do
   end do

    ! compute mtemp = A_0^n vbar^n
    call diffusive_m_fluxdiv(mla,mtemp,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^n vbar^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(mtemp(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)
       end do
    end do

    ! reset rho_update for all scalars to zero
    ! then, set rho_update to F^{*,n+1} = div(rho*chi grad c)^{*,n+1} + div(Psi^n)
    ! it is used in Step 5 below
    do n=1,nlevs
       call multifab_setval_c(rho_update(n),0.d0,1,nspecies,all=.true.)
       ! add fluxes
       call multifab_plus_plus_c(rho_update(n),1, diff_mass_fluxdiv(n),1,nspecies)
       if (variance_coef_mass .ne. 0.d0) then
          call multifab_plus_plus_c(rho_update(n),1,stoch_mass_fluxdiv(n),1,nspecies)
       end if
    end do

    ! compute div vbar^n
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div vbar^n to gmres_rhs_p
    ! now gmres_rhs_p = div vbar^n - S^{*,n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,1.d0/2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
       call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    gmres_abs_tol_in = gmres_abs_tol ! Save this 

    ! This relies entirely on relative tolerance and can fail if the rhs is roundoff error only:
    ! gmres_abs_tol = 0.d0 ! It is better to set gmres_abs_tol in namelist to a sensible value

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha,norm_pre_rhs)

    ! for the corrector gmres solve we want the stopping criteria based on the
    ! norm of the preconditioned rhs from the predictor gmres solve.  otherwise
    ! for cases where du in the corrector should be small the gmres stalls
    gmres_abs_tol = max(gmres_abs_tol_in, norm_pre_rhs*gmres_rel_tol)

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! compute v^{*,n+1} = v^n + dumac
    ! compute pi^{*,n+1}= pi^n + dpi
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_plus_plus_c(pi(n),1,dpi(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pi(n))
       call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! convert v^{*,n+1} to rho^{*,n+1}v^{*,n+1} in valid and ghost region
    ! now mnew has properly filled ghost cells
    call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 5 - Trapezoidal Scalar Corrector
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! rho_update already contains D^{*,n+1} + St^{*,n+1} for rho from above
    ! add A^{*,n+1} for rho to rho_update
    if (advection_type .ge. 1) then

       do n=1,nlevs
          ! bds force currently contains D^n + St^n
          ! add D^{*,n+1} + St^{*,n+1} and then multiply by 1/2
          call multifab_plus_plus_c(bds_force(n),1,rho_update(n),1,nspecies,0)
          call multifab_mult_mult_s_c(bds_force(n),1,0.5d0,nspecies,0)
          call multifab_fill_boundary(bds_force(n))
          do i=1,dm
             call multifab_plus_plus_c(umac_tmp(n,i),1,umac(n,i),1,1,1)
             call multifab_mult_mult_s_c(umac_tmp(n,i),1,0.5d0,1,1)
          end do
          call setval(rho_update(n),0.d0,all=.true.)
       end do

       if (advection_type .eq. 1 .or. advection_type .eq. 2) then

          call bds(mla,umac_tmp,rho_tmp,rho_update,bds_force,rho_fc,rho_nd_old,dx,dt,1,nspecies, &
                   c_bc_comp,the_bc_tower,proj_type_in=2)

          do n=1,nlevs
             call multifab_destroy(rho_nd_old(n))
             call multifab_destroy(rho_tmp(n))
          end do

       else if (advection_type .eq. 3 .or. advection_type .eq. 4) then
          call bds_quad(mla,umac_tmp,rho_old,rho_update,bds_force,rho_fc,dx,dt,1,nspecies, &
                        c_bc_comp,the_bc_tower,proj_type_in=2)
       end if    

       ! snew = s^n + dt * A^{n+1/2} + (dt/2) * (D^n + D^{n+1,*} + S^n + S^{n+1,*})
       do n=1,nlevs
          call multifab_copy_c(rho_new(n),1,rho_old(n),1,nspecies,0)
          call multifab_mult_mult_s_c(rho_update(n),1,dt,nspecies,0)
          call multifab_mult_mult_s_c(bds_force(n),1,dt,nspecies,0)
          call multifab_plus_plus_c(rho_new(n),1,rho_update(n),1,nspecies,0)
          call multifab_plus_plus_c(rho_new(n),1,bds_force(n),1,nspecies,0)
       end do

    else

       ! compute A^{*,n+1} for scalars
       call mk_advective_s_fluxdiv(mla,umac,rho_fc,rho_update,dx,1,nspecies)

       ! snew = s^{n+1} 
       !      = (1/2)*s^n + (1/2)*s^{*,n+1} + (dt/2)*(A^{*,n+1} + D^{*,n+1} + St^{*,n+1})
       do n=1,nlevs
          call multifab_plus_plus_c(rho_new(n),1,rho_old(n),1,nspecies,0)
          call multifab_mult_mult_s_c(rho_new(n),1,0.5d0,nspecies,0)
          call multifab_mult_mult_s_c(rho_update(n),1,dt/2.d0,nspecies,0)
          call multifab_plus_plus_c(rho_new(n),1,rho_update(n),1,nspecies,0)
       end do

    end if

    ! compute rhotot from rho in VALID REGION
    call compute_rhotot(mla,rho_new,rhotot_new)

    ! rho to conc - NO GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.true.)
    call fill_c_ghost_cells(mla,conc,dx,the_bc_tower)

    do n=1,nlevs
       call fill_rho_ghost_cells(conc(n),rhotot_new(n),the_bc_tower%bc_tower_array(n))
    end do

    ! conc to rho - INCLUDING GHOST CELLS
    call convert_rhoc_to_c(mla,rho_new,rhotot_new,conc,.false.)

    ! average rho_new and rhotot_new to faces
    call average_cc_to_face(nlevs,   rho_new,   rho_fc,1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array)
    call average_cc_to_face(nlevs,rhotot_new,rhotot_fc,1,    scal_bc_comp,       1,the_bc_tower%bc_tower_array)

    ! compute (eta,kappa)^{n+1}
    call compute_eta_kappa(mla,eta,eta_ed,kappa,rho_new,rhotot_new,Temp,dx, &
                           the_bc_tower%bc_tower_array)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Step 6 - Calculate Diffusive and Stochastic Fluxes
    ! Step 7 - Corrector Crank-Nicolson Step
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! build up rhs_v for gmres solve: first set gmres_rhs_v to mold/dt
    do n=1,nlevs
       do i=1,dm
          call multifab_copy_c(gmres_rhs_v(n,i),1,mold(n,i),1,1,0)
          call multifab_mult_mult_s_c(gmres_rhs_v(n,i),1,1.d0/dt,1,0)
       end do
    end do

    ! compute grad pi^{*,n+1}
    call compute_grad(mla,pi,gradpi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)

    ! barodiffusion uses predicted grad(pi)
    if (barodiffusion_type .eq. 2) then
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(gradp_baro(n,i),1,gradpi(n,i),1,1,0)
          end do
       end do
    else if (barodiffusion_type .eq. 3) then
       ! compute p0 from rho0*g
       call compute_HSE_pres(mla,rhotot_new,p_baro,dx,the_bc_tower)
       call compute_grad(mla,p_baro,gradp_baro,dx,1,pres_bc_comp,1,1, &
                         the_bc_tower%bc_tower_array)
    end if

    ! subtract grad pi^{*,n+1} from gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_sub_sub_c(gmres_rhs_v(n,i),1,gradpi(n,i),1,1,0)
       end do
    end do

    ! m_a_fluxdiv already contains A^n for momentum
    ! add A^{*,n+1} = -rho^{*,n+1} v^{*,n+1} v^{*,n+1} for momentum to m_a_fluxdiv
    call mk_advective_m_fluxdiv(mla,umac,mtemp,m_a_fluxdiv,dx, &
                                the_bc_tower%bc_tower_array)

    ! add (1/2) m_a_fluxdiv to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_a_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_a_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! m_d_fluxdiv already contains (1/2) A_0^n v^n
    ! add (1/2) A_0^n v^n to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute div(Sigma^n') by incrementing existing stochastic flux and dividing by 2
    if (variance_coef_mom .ne. 0.d0) then
       call stochastic_m_fluxdiv(mla,the_bc_tower%bc_tower_array,m_s_fluxdiv,eta,eta_ed, &
                                 Temp,Temp_ed,dx,dt,weights)
       do n=1,nlevs
          do i=1,dm
             call multifab_mult_mult_s_c(m_s_fluxdiv(n,i),1,0.5d0,1,0)
          end do
       end do
    end if

    ! add div(Sigma^n') to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_s_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! add gravity term
    if (any(grav(1:dm) .ne. 0.d0)) then
       call mk_grav_force(mla,gmres_rhs_v,rhotot_fc_old,rhotot_fc,the_bc_tower)
    end if

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,time+dt, &
                                   the_bc_tower%bc_tower_array)

    ! fill the stochastic multifabs with a new set of random numbers
    call fill_m_stochastic(mla)
    call fill_mass_stochastic(mla,the_bc_tower%bc_tower_array)

    ! compute diffusive and stochastic mass fluxes
    ! this computes "-F" so we later multiply by -1
    call compute_mass_fluxdiv(mla,rho_new,gradp_baro, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                      Temp,flux_total,dt,time,dx,weights, &
                                      the_bc_tower)

    do n=1,nlevs
       call multifab_mult_mult_s_c(diff_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       if (variance_coef_mass .ne. 0) then
          call multifab_mult_mult_s_c(stoch_mass_fluxdiv(n),1,-1.d0,nspecies,0)
       end if
       do i=1,dm
          call multifab_mult_mult_s_c(flux_total(n,i),1,-1.d0,nspecies,0)
       end do
    end do

    ! set the Dirichlet velocity value on reservoir faces
    call reservoir_bc_fill(mla,flux_total,vel_bc_n,the_bc_tower%bc_tower_array)

    ! compute gmres_rhs_p
    do n=1,nlevs
       call setval(gmres_rhs_p(n),0.d0,all=.true.)
       do i=1,nspecies
          call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call saxpy(gmres_rhs_p(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
       end do
    end do

    ! modify umac to respect the boundary conditions we want after the next gmres solve
    ! thus when we add A_0^{n+1} vbar^{*,n+1} to gmres_rhs_v and add div vbar^{*,n+1} to gmres_rhs_p we
    ! are automatically putting the system in delta form WITH homogeneous boundary conditions
    do n=1,nlevs
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

   ! compute mtemp = rho^{n+1} * vbar^{*,n+1}
   call convert_m_to_umac(mla,rhotot_fc,mtemp,umac,.false.)

   do n=1,nlevs
      do i=1,dm

         ! multiply mtemp by 1/dt
         call multifab_mult_mult_s_c(mtemp(n,i),1,1.d0/dt,1,0)

         ! subtract rho^{n+1} * vbar^{*,n+1} / dt from gmres_rhs_v
         call multifab_sub_sub_c(gmres_rhs_v(n,i),1,mtemp(n,i),1,1,0)

      end do
   end do

    ! set m_d_fluxdiv = A_0^{n+1} vbar^{n+1,*}
    do n=1,nlevs
       do i=1,dm
          call setval(m_d_fluxdiv(n,i),0.d0,all=.true.)
       end do
    end do
    call diffusive_m_fluxdiv(mla,m_d_fluxdiv,umac,eta,eta_ed,kappa,dx, &
                             the_bc_tower%bc_tower_array)

    ! add (1/2) A_0^{n+1} vbar^{n+1,*} to gmres_rhs_v
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_s_c(m_d_fluxdiv(n,i),1,0.5d0,1,0)
          call multifab_plus_plus_c(gmres_rhs_v(n,i),1,m_d_fluxdiv(n,i),1,1,0)
       end do
    end do

    ! compute div(vbar^{n+1,*})
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div(vbar^{n+1,*}) to gmres_rhs_p
    ! now gmres_rhs_p = div(vbar^{n+1,*}) - S^{n+1}
    ! the sign convention is correct since we solve -div(delta v) = gmres_rhs_p
    do n=1,nlevs
       call multifab_plus_plus_c(gmres_rhs_p(n),1,divu(n),1,1,0)
    end do

    ! multiply eta and kappa by 1/2 to put in proper form for gmres solve
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,1.d0/2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,1.d0/2.d0,1,eta(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,1.d0/2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! set the initial guess to zero
    do n=1,nlevs
       do i=1,dm
          call multifab_setval(dumac(n,i),0.d0,all=.true.)
       end do
          call multifab_setval(dpi(n),0.d0,all=.true.)
    end do

    do n=1,nlevs
       call zero_edgeval_physical(gmres_rhs_v(n,:),1,1,the_bc_tower%bc_tower_array(n))
    end do

    ! call gmres to compute delta v and delta pi
    call gmres(mla,the_bc_tower,dx,gmres_rhs_v,gmres_rhs_p,dumac,dpi,rhotot_fc, &
               eta,eta_ed,kappa,theta_alpha)
                              
    gmres_abs_tol = gmres_abs_tol_in ! Restore the desired tolerance   

    ! restore eta and kappa
    do n=1,nlevs
       call multifab_mult_mult_s_c(eta(n)  ,1,2.d0,1,eta(n)%ng)
       call multifab_mult_mult_s_c(kappa(n),1,2.d0,1,kappa(n)%ng)
       do i=1,size(eta_ed,dim=2)
          call multifab_mult_mult_s_c(eta_ed(n,i),1,2.d0,1,eta_ed(n,i)%ng)
       end do
    end do

    ! compute v^{n+1} = v^{n+1,*} + dumac
    ! compute pi^{n+1} = pi^{n+1,*} + dpi
    do n=1,nlevs
       do i=1,dm
          call multifab_plus_plus_c(umac(n,i),1,dumac(n,i),1,1,0)
       end do
       call multifab_plus_plus_c(pi(n),1,dpi(n),1,1,0)
    end do

    do n=1,nlevs
       ! presure ghost cells
       call multifab_fill_boundary(pi(n))
       call multifab_physbc(pi(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
       do i=1,dm
          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! set transverse velocity behind physical boundaries
          call multifab_physbc_macvel(umac(n,i),vel_bc_comp+i-1, &
                                      the_bc_tower%bc_tower_array(n), &
                                      dx(n,:),vel_bc_t(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! End Time-Advancement
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(rho_update(n))
       call multifab_destroy(bds_force(n))
       call multifab_destroy(gmres_rhs_p(n))
       call multifab_destroy(dpi(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(conc(n))
       call multifab_destroy(p_baro(n))
       do i=1,dm
          call multifab_destroy(mold(n,i))
          call multifab_destroy(mtemp(n,i))
          call multifab_destroy(m_a_fluxdiv(n,i))
          call multifab_destroy(m_d_fluxdiv(n,i))
          call multifab_destroy(m_s_fluxdiv(n,i))
          call multifab_destroy(gmres_rhs_v(n,i))
          call multifab_destroy(dumac(n,i))
          call multifab_destroy(umac_tmp(n,i))
          call multifab_destroy(rhotot_fc_old(n,i))
          call multifab_destroy(gradpi(n,i))
          call multifab_destroy(rho_fc(n,i))
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(flux_total(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine advance_timestep_inertial

  subroutine build_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs
    logical :: nodal_temp(3)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep_inertial/build_bc_multifabs")

    dm = mla%dim
    nlevs = mla%nlevel

    allocate(vel_bc_n(nlevs,dm))
    if (dm .eq. 2) then
       allocate(vel_bc_t(nlevs,2))
    else if (dm .eq. 3) then
       allocate(vel_bc_t(nlevs,6))
    end if

    do n=1,nlevs
       ! boundary conditions
       do i=1,dm
          call multifab_build_edge(vel_bc_n(n,i),mla%la(n),1,0,i)
       end do
       if (dm .eq. 2) then
          ! y-velocity bc on x-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,1),mla%la(n),1,0)
          ! x-velocity bc on y-faces (nodal)
          call multifab_build_nodal(vel_bc_t(n,2),mla%la(n),1,0)
       else
          ! y-velocity bc on x-faces (nodal in y and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,1),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on x-faces (nodal in z and x)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,2),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on y-faces (nodal in x and y)
          nodal_temp(1) = .true.
          nodal_temp(2) = .true.
          nodal_temp(3) = .false.
          call multifab_build(vel_bc_t(n,3),mla%la(n),1,0,nodal_temp)
          ! z-velocity bc on y-faces (nodal in z and y)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,4),mla%la(n),1,0,nodal_temp)
          ! x-velocity bc on z-faces (nodal in x and z)
          nodal_temp(1) = .true.
          nodal_temp(2) = .false.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,5),mla%la(n),1,0,nodal_temp)
          ! y-velocity bc on z-faces (nodal in y and z)
          nodal_temp(1) = .false.
          nodal_temp(2) = .true.
          nodal_temp(3) = .true.
          call multifab_build(vel_bc_t(n,6),mla%la(n),1,0,nodal_temp)
       end if

       do i=1,dm
          call multifab_setval(vel_bc_n(n,i),0.d0,all=.true.)
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_setval(vel_bc_t(n,i),0.d0,all=.true.)
       end do

    end do

    call destroy(bpt)

  end subroutine build_bc_multifabs

  subroutine destroy_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt, "advance_timestep_inertial/destroy_bc_multifabs")

    dm = mla%dim
    nlevs = mla%nlevel

    do n=1,nlevs
       do i=1,dm          
          call multifab_destroy(vel_bc_n(n,i))
       end do
       do i=1,size(vel_bc_t,dim=2)
          call multifab_destroy(vel_bc_t(n,i))
       end do
    end do

    deallocate(vel_bc_n,vel_bc_t)

    call destroy(bpt)

  end subroutine destroy_bc_multifabs

end module advance_timestep_inertial_module
