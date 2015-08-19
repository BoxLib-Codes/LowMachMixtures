module apply_precon_module

  use ml_layout_module
  use multifab_module
  use define_bc_module
  use probin_gmres_module , only: precon_type, visc_schur_approx
  use probin_common_module, only: visc_type
  use stag_mg_solver_module
  use stag_applyop_module
  use macproject_module
  use div_and_grad_module
  use cc_applyop_module
  use convert_stag_module
  use norm_inner_product_module
  use multifab_physbc_module
  use bc_module
  use inverse_diag_lap_module
  use stag_mg_layout_module

  implicit none

  private

  public :: apply_precon

contains

  ! This computes x = M^{-1} b using the approach in ./doc/PreconditionerNotes.tex
  subroutine apply_precon(mla,b_u,b_p,x_u,x_p,alpha_fc,beta,beta_ed, &
                          gamma,theta_alpha,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: b_u(:,:)
    type(multifab) , intent(in   ) :: b_p(:)
    type(multifab) , intent(inout) :: x_u(:,:)
    type(multifab) , intent(inout) :: x_p(:)
    type(multifab) , intent(in   ) :: alpha_fc(:,:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: beta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma(:)
    real(kind=dp_t), intent(in   ) :: theta_alpha
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs,i,lohi,dm
    real(kind=dp_t) :: mean_val_pres, mean_val_umac(mla%dim)

    type(multifab) ::         phi(mla%nlevel)
    type(multifab) ::     mac_rhs(mla%nlevel)
    type(multifab) ::    zero_fab(mla%nlevel)
    type(multifab) ::     x_p_tmp(mla%nlevel)
    type(multifab) :: alphainv_fc(mla%nlevel,mla%dim)
    type(multifab) ::     b_u_tmp(mla%nlevel,mla%dim)
    type(multifab) ::  one_fab_fc(mla%nlevel,mla%dim)
    type(multifab) :: zero_fab_fc(mla%nlevel,mla%dim)

    logical :: no_wall_is_no_slip

    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "apply_precon")

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(phi(n)     ,mla%la(n),1,1)
       call multifab_build(mac_rhs(n) ,mla%la(n),1,0)
       call multifab_build(zero_fab(n),mla%la(n),1,0)
       call setval(zero_fab(n),0.d0,all=.true.)
       call multifab_build(x_p_tmp(n),mla%la(n),1,1)
       do i=1,dm
          call multifab_build_edge(alphainv_fc(n,i),mla%la(n),1,0,i)
          call setval(alphainv_fc(n,i),1.d0,all=.true.)
          call multifab_div_div_c(alphainv_fc(n,i),1,alpha_fc(n,i),1,1,0)
          call multifab_build_edge(b_u_tmp(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(one_fab_fc(n,i),mla%la(n),1,0,i)
          call setval(one_fab_fc(n,i),1.d0,all=.true.)
          call multifab_build_edge(zero_fab_fc(n,i),mla%la(n),1,0,i)
          call setval(zero_fab_fc(n,i),0.d0,all=.true.)
       end do
    end do

    ! set the initial guess for Phi in the Poisson solve to 0
    ! set x_u = 0 as initial guess
    do n=1,nlevs
       call setval(phi(n),0.d0,all=.true.)
       do i=1,dm
          call setval(x_u(n,i),0.d0,all=.true.)
       end do
    end do

    ! 1 = projection preconditioner
    ! 2 = lower triangular preconditioner
    ! 3 = upper triangular preconditioner
    ! 4 = block diagonal preconditioner
    ! 5 = Uzawa-type approximation (see paper)
    ! 6 = upper triangular + viscosity-based BFBt Schur complement (from Georg Stadler)

    select case (abs(precon_type))

    case(1) ! projection preconditioner

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! STEP 1: Solve for an intermediate state, x_u^star, using an implicit viscous solve
       !         x_u^star = A^{-1} b_u
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! x_u^star = A^{-1} b_u
        call stag_mg_solver(mla,la_mg_normal,alpha_fc,beta,beta_ed,gamma,theta_alpha,x_u,b_u,dx,the_bc_tower)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Construct RHS for pressure Poisson problem
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! set mac_rhs = D(x_u^star)
        call compute_div(mla,x_u,mac_rhs,dx,1,1,1)

        ! add b_p to mac_rhs
        do n=1,nlevs
          call multifab_plus_plus_c(mac_rhs(n),1,b_p(n),1,1,0)
        end do

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Compute x_u
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! use multigrid to solve for Phi
        ! x_u^star is only passed in to get a norm for absolute residual criteria
        call macproject(mla,phi,x_u,alphainv_fc,mac_rhs,dx,the_bc_tower)

        ! x_u = x_u^star - (alpha I)^-1 grad Phi
        call subtract_weighted_gradp(mla,x_u,alphainv_fc,phi,dx,the_bc_tower)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 4: Compute x_p by applying the Schur complement approximation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (visc_schur_approx .eq. 0) then
           
           ! if precon_type = +1, or theta_alpha=0 then x_p = theta_alpha*Phi - c*beta*(mac_rhs)
           ! if precon_type = -1                   then x_p = theta_alpha*Phi - c*beta*L_alpha Phi
           
           if ((precon_type .eq. 1) .or. (theta_alpha .eq. 0.d0)) then
              ! first set x_p = -mac_rhs 
              do n=1,nlevs
                 call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)
                 call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
              end do
           else
              ! first set x_p = -L_alpha Phi
              call cc_applyop(mla,x_p,phi,zero_fab,alphainv_fc,dx, &
                              the_bc_tower,pres_bc_comp,stencil_order_in=2)
           end if

           do n=1,nlevs

              if ( (abs(visc_type) .eq. 1) .or. (abs(visc_type) .eq. 2) ) then  
                 ! multiply x_p by beta; x_p = -beta L_alpha Phi
                 call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)

                 if (abs(visc_type) .eq. 2) then
                    ! multiply by c=2; x_p = -2*beta L_alpha Phi
                    call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
                 end if

              else if (abs(visc_type) .eq. 3) then

                 ! multiply x_p by gamma, use mac_rhs a temparary to save x_p 
                 call multifab_copy_c(mac_rhs(n),1,x_p(n),1,1,0)
                 call multifab_mult_mult_c(mac_rhs(n),1,gamma(n),1,1,0)
                 ! multiply x_p by beta; x_p = -beta L_alpha Phi
                 call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
                 ! multiply by c=4/3; x_p = -(4/3) beta L_alpha Phi
                 call multifab_mult_mult_s_c(x_p(n),1,4.d0/3.d0,1,0)
                 ! x_p = -(4/3) beta L_alpha Phi - gamma L_alpha Phi
                 call multifab_plus_plus_c(x_p(n),1,mac_rhs(n),1,1,0)

              end if

              ! multiply Phi by theta_alpha
              call multifab_mult_mult_s_c(phi(n),1,theta_alpha,1,0)

              ! add theta_alpha*Phi to x_p
              call multifab_plus_plus_c(x_p(n),1,phi(n),1,1,0)

           end do

        else

           call visc_schur_complement(mla,mac_rhs,x_p,x_u,beta,beta_ed, &
                                      theta_alpha,dx,the_bc_tower,phi)
           
        end if

      case(2,5)  
        ! lower triangular, precon_type=-2 means using negative sign for Schur complement app

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve for x_u using an implicit viscous term
        !         A x_u = b_u
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! x_u = A^{-1} b_u
        call stag_mg_solver(mla,la_mg_normal,alpha_fc,beta,beta_ed,gamma,theta_alpha,x_u,b_u,dx,the_bc_tower)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Solve a pressure Poisson problem for Phi
        !         L_alpha Phi = D x_u + b_p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! set mac_rhs = D(x_u)
        call compute_div(mla,x_u,mac_rhs,dx,1,1,1)

        ! add b_p to mac_rhs
        do n=1,nlevs
          call multifab_plus_plus_c(mac_rhs(n),1,b_p(n),1,1,0)
        end do

        if (abs(theta_alpha) .gt. 0.d0) then

           ! solves L_alpha Phi = mac_rhs
           ! x_u is only passed in to get a norm for absolute residual criteria
           call macproject(mla,phi,x_u,alphainv_fc,mac_rhs,dx,the_bc_tower)

        end if 
    
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Compute x_p by applying the Schur complement approximation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (visc_schur_approx .eq. 0) then
         
           call simple_schur_complement()

        else

           call visc_schur_complement(mla,mac_rhs,x_p,x_u,beta,beta_ed, &
                                      theta_alpha,dx,the_bc_tower,phi)

        end if

        do n=1,nlevs 
           if (precon_type .eq. -2 .or. precon_type .eq. -5) then
              ! multiply x_p by -1, if precon_type=-2
              call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
           end if
        end do   

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 5: Additional steps for P_5
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (abs(precon_type) .eq. 5) then

           ! compute = A^(-1)*(b_u-grad(x_p)) 

           ! we need gradients of x_p, and x_p doesn't necessarily have 
           ! a ghost cell so we create a temporary
           do n=1,nlevs
              call multifab_copy_c(x_p_tmp(n),1,x_p(n),1,1,0)
              call multifab_fill_boundary(x_p_tmp(n))
              call multifab_physbc(x_p_tmp(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                                   dx_in=dx(n,:))
           end do

           ! contstruct b_u-grad(x_p)
           do n=1,nlevs
              do i=1,dm
                 call multifab_copy_c(b_u_tmp(n,i),1,b_u(n,i),1,1,0)
              end do
           end do
           call subtract_weighted_gradp(mla,b_u_tmp,one_fab_fc,x_p_tmp,dx,the_bc_tower)

           ! compute = A^(-1)*(b_u-grad(x_p)) 
           call stag_mg_solver(mla,la_mg_normal,alpha_fc,beta,beta_ed,gamma,theta_alpha,x_u,b_u_tmp, &
                               dx,the_bc_tower)

        end if
        
      case(3)  
        ! upper triangular, precon_type=-3 means using negative sign for Schur complement app

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve a pressure Poisson problem for Phi
        !         L_alpha Phi = b_p
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do n=1,nlevs
          ! copy b_p to mac_rhs, mac_rhs is now the rhs for later use 
          call multifab_copy_c(mac_rhs(n),1,b_p(n),1,1,0)  
          if (precon_type .eq. -3) then         ! rhs times -1   
            call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
          end if 
        end do
  
        if (abs(theta_alpha) .gt. 0) then

          ! solves L_alpha Phi = mac_rhs
          ! x_u^star is only passed in to get a norm for absolute residual criteria
          call macproject(mla,phi,x_u,alphainv_fc,mac_rhs,dx,the_bc_tower)

        end if
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Compute x_p by applying the Schur complement approximation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (visc_schur_approx .eq. 0) then
           
           call simple_schur_complement()

        else

           call visc_schur_complement(mla,mac_rhs,x_p,x_u,beta,beta_ed, &
                                      theta_alpha,dx,the_bc_tower,phi)

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 3: Compute the RHS for the viscous solve
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! we need gradients of x_p, and x_p doesn't necessarily have 
        ! a ghost cell so we create a temporary
        do n=1,nlevs
          call multifab_copy_c(x_p_tmp(n),1,x_p(n),1,1,0)
          call multifab_fill_boundary(x_p_tmp(n))
          call multifab_physbc(x_p_tmp(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                               dx_in=dx(n,:))
        end do

        ! contstruct b_u-grad(x_p)
        do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(b_u_tmp(n,i),1,b_u(n,i),1,1,0)
          end do
        end do
        call subtract_weighted_gradp(mla,b_u_tmp,one_fab_fc,x_p_tmp,dx,the_bc_tower)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 4: Solve for x_u using an implicit viscous term
        !         A x_u = (b_u-G x_p)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! compute = A^(-1)*(b_u-grad(x_p)) 
        call stag_mg_solver(mla,la_mg_normal,alpha_fc,beta,beta_ed,gamma,theta_alpha,x_u,b_u_tmp, &
                            dx,the_bc_tower)

      case(4) 
        ! block diagonal

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 1: Solve for x_u using an implicit viscous term
        !         A x_u = b_u
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! x_u = A^{-1} b_u
        call stag_mg_solver(mla,la_mg_normal,alpha_fc,beta,beta_ed,gamma,theta_alpha,x_u,b_u,dx,the_bc_tower)

        if (abs(theta_alpha) .gt. 0) then  

          do n=1,nlevs
            ! copy b_p to mac_rhs, mac_rhs is now the rhs for later use 
            call multifab_copy_c(mac_rhs(n),1,b_p(n),1,1,0)  
          end do

          ! solves L_alpha Phi = mac_rhs
          ! x_u^star is only passed in to get a norm for absolute residual criteria
          call macproject(mla,phi,x_u,alphainv_fc,mac_rhs,dx,the_bc_tower)

        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! STEP 2: Compute x_p by applying the Schur complement approximation
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (visc_schur_approx .eq. 0) then
           
           do n=1,nlevs

              if (abs(theta_alpha) .gt. 0) then
                 call multifab_mult_mult_s_c(phi(n),1,theta_alpha,1,0)
              end if
              ! x_p = theta_alpha*Phi-beta*bp
              call multifab_copy_c(x_p(n),1,b_p(n),1,1,0)
              call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)
              call multifab_sub_sub_c(x_p(n),1,phi(n),1,1,0) 

              if ((abs(visc_type) .eq. 1) .and. (precon_type .eq. 4)) then
                 ! multiply by -c=-1 for |viscous_type| = 1
                 call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)
              else if ((abs(visc_type) .eq. 2) .and. (precon_type .eq. 4)) then
                 ! multiply by -c=-2 for |viscous_type| = 2
                 call multifab_mult_mult_s_c(x_p(n),1,-2.d0,1,0)
              else if ((abs(visc_type) .eq. 2) .and. (precon_type .eq. -4)) then
                 ! multiply by c=2 for precon_type = -4
                 call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
              else if ((abs(visc_type) .eq. 3) .and. (precon_type .eq. 4)) then
                 ! multiply by -c=-4/3 for |viscous_type| = 3
                 call multifab_mult_mult_s_c(x_p(n),1,-4.d0/3.d0,1,0)
                 ! gamma part: the sign is same as beta, use x_p_tmp as an immediate variable  
                 call multifab_copy_c(x_p_tmp(n),1,b_p(n),1,1,0)
                 call multifab_mult_mult_c(x_p_tmp(n),1,gamma(n),1,1,0)
                 call multifab_mult_mult_s_c(x_p_tmp(n),1,-1.d0,1,0)
                 call multifab_plus_plus_c(x_p(n),1,x_p_tmp(n),1,1,0)
              else if ((abs(visc_type) .eq. 3) .and. (precon_type .eq. -4)) then
                 ! multiply by c=4/3 for precon_type = -4
                 call multifab_mult_mult_s_c(x_p(n),1,4.d0/3.d0,1,0)
                 ! gamma part: the sign is same as beta, use x_p_tmp as an immediate variable   
                 call multifab_copy_c(x_p_tmp(n),1,b_p(n),1,1,0)
                 call multifab_mult_mult_c(x_p_tmp(n),1,gamma(n),1,1,0)
                 call multifab_plus_plus_c(x_p(n),1,x_p_tmp(n),1,1,0)             
              end if

           end do

        else

           call visc_schur_complement(mla,mac_rhs,x_p,x_u,beta,beta_ed, &
                                      theta_alpha,dx,the_bc_tower,phi)

        end if

     case default
        call bl_error('apply_precon.f90: unsupported precon_type')
     end select

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 4: Handle null-space issues in MG solvers
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! subtract off mean value: Single level only! No need for ghost cells
    call sum_umac_press(mla,x_p,x_u,mean_val_pres,mean_val_umac)
    
    ! The pressure Poisson problem is always singular:
    call multifab_sub_sub_s_c(x_p(1),1,mean_val_pres,1,0)
    
    ! The velocity problem is also singular under these cases
    if (theta_alpha .eq. 0.d0) then

       no_wall_is_no_slip = .true.
       do i=1,dm
          do lohi=1,2
             if (the_bc_tower%domain_bc(i,lohi) .ge. NO_SLIP_START .and. &
                 the_bc_tower%domain_bc(i,lohi) .le. NO_SLIP_END) then
                 no_wall_is_no_slip = .false.
              end if
          end do
       end do

       if (no_wall_is_no_slip) then
          do i=1,dm
             if (mla%pmask(i)) then
                call multifab_sub_sub_s_c(x_u(1,i),1,mean_val_umac(i),1,0)
             end if
          end do
       end if

    end if

    do n=1,nlevs
       call multifab_destroy(phi(n))
       call multifab_destroy(mac_rhs(n))
       call multifab_destroy(zero_fab(n))
       call multifab_destroy(x_p_tmp(n))
       do i=1,dm
          call multifab_destroy(alphainv_fc(n,i))
          call multifab_destroy(one_fab_fc(n,i))
          call multifab_destroy(zero_fab_fc(n,i))
          call multifab_destroy(b_u_tmp(n,i))
       end do
    end do

    call destroy(bpt)

  contains
  
    subroutine simple_schur_complement()
  
      type(bl_prof_timer), save :: bpt2
    
      call build(bpt2, "simple_schur_complement")

       ! x_p = theta_alpha*I *Phi - beta*mac_rhs 
       do n=1,nlevs
          ! beta part
          if ( (abs(visc_type) .eq. 1) .or. (abs(visc_type) .eq. 2) ) then 
             call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)

             ! multiply x_p by -1
             call multifab_mult_mult_s_c(x_p(n),1,-1.d0,1,0)

             ! multiply x_p by beta; x_p = -beta L_alpha Phi
             call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)

             if (abs(visc_type) .eq. 2) then
                ! multiply by c=2 for |viscous_type| = 2
                call multifab_mult_mult_s_c(x_p(n),1,2.d0,1,0)
             end if

          elseif (abs(visc_type) .eq. 3) then   
             ! beta part 
             call multifab_copy_c(x_p(n),1,mac_rhs(n),1,1,0)

             ! multiply x_p by -4/3
             call multifab_mult_mult_s_c(x_p(n),1,-4.d0/3.d0,1,0)

             ! multiply x_p by beta; x_p = -beta L_alpha Phi
             call multifab_mult_mult_c(x_p(n),1,beta(n),1,1,0)

             ! gamma part
             ! x_p = theta_alpha*I *Phi - 4/3*beta*mac_rhs - gamma*mac_rhs 
             call multifab_mult_mult_s_c(mac_rhs(n),1,-1.d0,1,0)
             call multifab_mult_mult_c(mac_rhs(n),1,gamma(n),1,1,0)

             ! x_p = x_p + mac_rhs 
             call multifab_plus_plus_c(x_p(n),1,mac_rhs(n),1,1,0)

          end if

          if (abs(theta_alpha) .gt. 0.d0) then         
             ! multiply Phi by theta_alpha
             call multifab_mult_mult_s_c(phi(n),1,theta_alpha,1,0)   

             ! add theta_alpha*Phi to x_p
             call multifab_plus_plus_c(x_p(n),1,phi(n),1,1,0)
          end if

       end do
    
       call destroy(bpt2)

      end subroutine simple_schur_complement
      
  end subroutine apply_precon

  subroutine visc_schur_complement(mla,rhs,x_out,x_u,beta,beta_ed, &
                                   theta_alpha,dx,the_bc_tower,phi)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rhs(:)
    type(multifab) , intent(inout) :: x_out(:)
    type(multifab) , intent(in   ) :: x_u(:,:)
    type(multifab) , intent(in   ) :: beta(:)
    type(multifab) , intent(in   ) :: beta_ed(:,:)
    real(kind=dp_t), intent(in   ) :: theta_alpha
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: phi(:)

    ! local
    integer:: i,dm,n,nlevs

    type(multifab) ::      phi_mu(mla%nlevel)
    type(multifab) ::     phi_rho(mla%nlevel)
    type(multifab) ::    zero_fab(mla%nlevel)
    type(multifab) ::        gphi(mla%nlevel,mla%dim)
    type(multifab) ::       Lgphi(mla%nlevel,mla%dim)
    type(multifab) ::    muinv_fc(mla%nlevel,mla%dim)
    type(multifab) :: zero_fab_fc(mla%nlevel,mla%dim)

    type(bl_prof_timer), save :: bpt
    
    call build(bpt, "visc_schur_complement")

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       call multifab_build(phi_mu(n),mla%la(n),1,1)
       call multifab_build(phi_rho(n),mla%la(n),1,1)
       call multifab_build(zero_fab(n),mla%la(n),1,0)
       call setval(zero_fab(n),0.d0,all=.true.)
       do i=1,dm
          call multifab_build_edge(gphi(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(Lgphi(n,i),mla%la(n),1,1,i)
          call multifab_build_edge(muinv_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(zero_fab_fc(n,i),mla%la(n),1,0,i)
          call setval(zero_fab_fc(n,i),0.d0,all=.true.)
       end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 1: Compute L_invmu^(-1)*rhs = (D*M^(-1)*G)^(-1)*rhs
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! compute coefficients on edges using M, the digonal of the viscous operator
    call inverse_diag_lap(mla,beta,beta_ed,muinv_fc)

    ! solve for a new Phi with inverse-viscosity weighted Poisson solve
    ! first reset initial guess for phi to zero
    do n=1,nlevs
       call multifab_setval(phi_mu(n),0.d0,all=.true.)
    end do

    ! x_u^star is only passed in to get a norm for absolute residual criteria
    call macproject(mla,phi_mu,x_u,muinv_fc,rhs,dx,the_bc_tower)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 2: Multiply result by D*M^(-1)*L_mu*M^(-1)*G
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! take gradient of new Phi
    call compute_grad(mla,phi_mu,gphi,dx,1,pres_bc_comp,1,1,the_bc_tower%bc_tower_array)
    
    ! multiply gradient by face-centered inverse diagonal coefficient
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(gphi(n,i),1,muinv_fc(n,i),1,1,0)
          call multifab_fill_boundary(gphi(n,i))
       end do
    end do

    ! apply the viscous operator
    call stag_applyop(mla,the_bc_tower,gphi,Lgphi,zero_fab_fc, &
                      beta,beta_ed,zero_fab,theta_alpha,dx)

    ! multiply result by face-centered inverse diagonal coefficient
    do n=1,nlevs
       do i=1,dm
          call multifab_mult_mult_c(Lgphi(n,i),1,muinv_fc(n,i),1,1,0)
          call multifab_fill_boundary(Lgphi(n,i))
       end do
    end do

    ! take divergence and store in rhs
    call compute_div(mla,Lgphi,rhs,dx,1,1,1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 3: Multiply result by L_invmu^(-1) = (D*M^(-1)*G)^(-1)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! solve an inverse-viscosity weighted Poisson solve
    ! first reset initial guess for phi to zero
    do n=1,nlevs
       call multifab_setval(phi_mu(n),0.d0,all=.true.)
    end do

    ! x_u^star is only passed in to get a norm for absolute residual criteria
    call macproject(mla,phi_mu,x_u,muinv_fc,rhs,dx,the_bc_tower)

    ! copy solution of Poisson equation into x_out
    do n=1,nlevs
       call multifab_copy_c(x_out(n),1,phi_mu(n),1,1,0)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! STEP 4: Add inertial correction if any
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (abs(theta_alpha) .gt. 0) then  

       ! we have already computed phi=Ltilde_rho^inv*rhs
       do n=1,nlevs
          call multifab_copy_c(phi_rho(n),1,phi(n),1,1,0)
       end do
    
       ! multiply phi_rho by theta_alpha
       do n=1,nlevs
          call multifab_mult_mult_s_c(phi_rho(n),1,theta_alpha,1,0)
       end do

       ! add theta_alpha*phi_rho to x_out
       do n=1,nlevs
          call multifab_plus_plus_c(x_out(n),1,phi_rho(n),1,1,0)
       end do
    
    end if   

    do n=1,nlevs
       call multifab_destroy(phi_mu(n))
       call multifab_destroy(phi_rho(n))
       call multifab_destroy(zero_fab(n))
       do i=1,dm
          call multifab_destroy(gphi(n,i))
          call multifab_destroy(Lgphi(n,i))
          call multifab_destroy(muinv_fc(n,i))
          call multifab_destroy(zero_fab_fc(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine visc_schur_complement

end module apply_precon_module
