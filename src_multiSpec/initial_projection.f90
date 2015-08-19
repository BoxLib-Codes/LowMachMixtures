module initial_projection_module

  use multifab_module
  use ml_layout_module
  use convert_stag_module
  use define_bc_module
  use macproject_module
  use div_and_grad_module
  use bc_module
  use multifab_physbc_stag_module
  use compute_mass_fluxdiv_module
  use reservoir_bc_fill_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: rhobar, variance_coef_mass, algorithm_type

  implicit none

  private

  public :: initial_projection

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

  subroutine initial_projection(mla,umac,rho,rhotot,gradp_baro, &
                                diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                Temp,eta,eta_ed,dt,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: umac(:,:)
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(inout) :: rhotot(:)
    type(multifab) , intent(in   ) :: gradp_baro(:,:)
    type(multifab) , intent(inout) :: diff_mass_fluxdiv(:)
    type(multifab) , intent(inout) :: stoch_mass_fluxdiv(:)
    type(multifab) , intent(in   ) :: Temp(:)
    type(multifab) , intent(in   ) :: eta(:)
    type(multifab) , intent(in   ) :: eta_ed(:,:)  ! nodal (2d); edge-centered (3d)
    real(kind=dp_t), intent(in   ) :: dt
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,dm,n,nlevs

    type(multifab) ::mac_rhs(mla%nlevel)
    type(multifab) :: divu(mla%nlevel)
    type(multifab) :: phi(mla%nlevel)
    type(multifab) :: rhotot_fc(mla%nlevel,mla%dim)
    type(multifab) :: rhototinv_fc(mla%nlevel,mla%dim)
    type(multifab) :: flux_total(mla%nlevel,mla%dim)

    real(kind=dp_t) :: weights(algorithm_type)
    
    type(bl_prof_timer), save :: bpt

    call build(bpt,"initial_projection")

    if (algorithm_type .eq. 1) then
       weights(1) = 1.d0
    else if (algorithm_type .eq. 2) then
       weights(1) = 1.d0
       weights(2) = 0.d0
    end if

    dm = mla%dim
    nlevs = mla%nlevel
 
    call build_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_build(mac_rhs(n),mla%la(n),1,0)
       call multifab_build(divu(n),mla%la(n),1,0)
       call multifab_build(phi(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(   rhotot_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(rhototinv_fc(n,i),mla%la(n),1       ,0,i)
          call multifab_build_edge(  flux_total(n,i),mla%la(n),nspecies,0,i)
       end do       
    end do

    do n=1,nlevs
       call setval(mac_rhs(n),0.d0)
       call setval(phi(n),0.d0)
    end do

    ! reset inhomogeneous bc condition to deal with reservoirs
    call set_inhomogeneous_vel_bcs(mla,vel_bc_n,vel_bc_t,eta_ed,dx,0.d0, &
                                   the_bc_tower%bc_tower_array)

    ! compute diff/stoch/baro mass fluxes
    ! this computes "-F" so we later multiply by -1
    call compute_mass_fluxdiv(mla,rho,gradp_baro, &
                                      diff_mass_fluxdiv,stoch_mass_fluxdiv, &
                                      Temp, &
                                      flux_total,dt,0.d0,dx,weights, &
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

    ! set mac_rhs to -S
    ! -S = -sum div (F_i / rhobar_i)
    do n=1,nlevs
       do i=1,nspecies
          call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i), diff_mass_fluxdiv(n),i,1)
          if (variance_coef_mass .ne. 0.d0) then
             call multifab_saxpy_3_cc(mac_rhs(n),1,-1.d0/rhobar(i),stoch_mass_fluxdiv(n),i,1)
          end if
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! build rhs = div(v^init) - S^0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    do n=1,nlevs
       do i=1,dm
!          ! to deal with reservoirs
!          ! set normal velocity on physical domain boundaries
          call multifab_physbc_domainvel(umac(n,i),vel_bc_comp+i-1, &
                                         the_bc_tower%bc_tower_array(n), &
                                         dx(n,:),vel_bc_n(n,:))
          ! fill periodic and interior ghost cells
          call multifab_fill_boundary(umac(n,i))
       end do
    end do

    ! set divu = div(v^init)
    call compute_div(mla,umac,divu,dx,1,1,1)

    ! add div(v^init) to mac_rhs
    ! now mac_rhs = div(v^init) - S
    do n=1,nlevs
       call multifab_plus_plus_c(mac_rhs(n),1,divu(n),1,1,0)
    end do

    ! average rhotot to faces
    call average_cc_to_face(nlevs,rhotot,rhotot_fc,1,scal_bc_comp,1,the_bc_tower%bc_tower_array)

    ! compute (1/rhotot) on faces
    do n=1,nlevs
       do i=1,dm
          call setval(rhototinv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(rhototinv_fc(n,i),1,rhotot_fc(n,i),1,1,0)
       end do
    end do

    ! solve div (1/rhotot) grad phi = div(v^init) - S^0
    ! solve to completion, i.e., use the 'full' solver
    call macproject(mla,phi,umac,rhototinv_fc,mac_rhs,dx,the_bc_tower,.true.)

    ! v^0 = v^init - (1/rho^0) grad phi
    call subtract_weighted_gradp(mla,umac,rhototinv_fc,phi,dx,the_bc_tower)

    ! fill ghost cells
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

    call destroy_bc_multifabs(mla)

    do n=1,nlevs
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(divu(n))
       call multifab_destroy(phi(n))
       do i=1,dm
          call multifab_destroy(rhotot_fc(n,i))
          call multifab_destroy(rhototinv_fc(n,i))
          call multifab_destroy(flux_total(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine initial_projection

  subroutine build_bc_multifabs(mla)

    type(ml_layout), intent(in   ) :: mla

    integer :: dm,i,n,nlevs
    logical :: nodal_temp(3)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"initial_projection/build_bc_multifabs")

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

    call build(bpt,"initial_projection/destroy_bc_multifabs")

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

end module initial_projection_module
