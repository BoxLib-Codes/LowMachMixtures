!**************************************************************************
! Compute the solution to a nonsymmetric system of linear equations of the
! form Ax=b by the restarted left-preconditioned Generalized Minimal Residual method (GMRES), 
! This code is based on Givens tranform, not Housholder transform, we refer to the Matlab code in 
! http://www.netlib.org/templates/matlab/gmres.m
!**************************************************************************
module gmres_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bl_constants_module
  use bc_module
  use bl_types
  use probin_gmres_module, only: gmres_max_inner, gmres_abs_tol, gmres_max_iter, &
       gmres_max_outer, gmres_min_iter, gmres_rel_tol, gmres_verbose, p_norm_weight, &
       scale_factor
  use vcycle_counter_module
  use apply_precon_module
  use apply_matrix_module
  use norm_inner_product_module

  implicit none 
   
  private
  public :: gmres
   
contains
  
  subroutine gmres(mla,the_bc_tower,dx,b_u,b_p,x_u,x_p, &
                   alpha_fc,beta,beta_ed,gamma,theta_alpha,norm_pre_rhs)

    type(ml_layout),intent(in   ) :: mla
    type(bc_tower), intent(in   ) :: the_bc_tower
    real(dp_t)    , intent(in   ) :: dx(:,:)
    type(multifab), intent(inout) :: b_u(:,:)  ! (nlevs,dm)
    type(multifab), intent(in   ) :: b_p(:)    ! (nlevs)
    type(multifab), intent(inout) :: x_u(:,:)  ! (nlevs,dm)
    type(multifab), intent(inout) :: x_p(:)    ! (nlevs)
    type(multifab), intent(in   ) :: alpha_fc(:,:)
    type(multifab), intent(inout) :: beta(:)
    type(multifab), intent(inout) :: beta_ed(:,:) ! nodal (2d), edge-centered (3d)
    type(multifab), intent(inout) :: gamma(:)
    real(dp_t)    , intent(inout) :: theta_alpha
    real(dp_t), optional, intent(out) :: norm_pre_rhs

    ! Local
    type(multifab) ::   r_u(mla%nlevel,mla%dim)
    type(multifab) ::   w_u(mla%nlevel,mla%dim)
    type(multifab) :: tmp_u(mla%nlevel,mla%dim)
    type(multifab) ::   r_p(mla%nlevel)
    type(multifab) ::   w_p(mla%nlevel)
    type(multifab) :: tmp_p(mla%nlevel)

    real(kind=dp_t) :: cs(gmres_max_inner)
    real(kind=dp_t) :: sn(gmres_max_inner)
    real(kind=dp_t) ::  y(gmres_max_inner)
    real(kind=dp_t) ::  s(gmres_max_inner+1)
    real(kind=dp_t) ::  H(gmres_max_inner+1,gmres_max_inner) ! for storing Hessenberg mat

    type(bl_prof_timer), save :: bpt

    ! Krylov space vectors
    type(multifab) :: V_u(mla%nlevel,mla%dim)
    type(multifab) :: V_p(mla%nlevel)

    integer :: n, d, i, k, dm, nlevs, iter, total_iter, i_copy  ! for looping iteration
    integer :: vcycle_counter_temp ! Local copy
    real(kind=dp_t) :: norm_b,norm_pre_b,norm_resid,rel_resid, norm_init_resid
    real(kind=dp_t) :: norm_resid_Stokes, norm_init_Stokes
    real(kind=dp_t) :: norm_u, norm_p, norm_init_resid_u, norm_init_resid_p
    real(kind=dp_t) :: norm_init_u, norm_init_p, norm_u_noprecon
    real(kind=dp_t) :: norm_p_noprecon, norm_resid_est
    real(kind=dp_t) :: inner_prod_vel(mla%dim), inner_prod_pres

    call build(bpt, "gmres")

    if (gmres_verbose .ge. 3) then 
      if (parallel_IOProcessor()) then 
        open (unit = 204, file = "gmres_inner_resid.dat", form="formatted", status="unknown", action = "write")
        write ( 204, '(A)' ) "# 1=#V, 2=|Pr_est|/|Pr_est0|"

        open (unit = 205, file = "gmres_outer_resid.dat", form="formatted", status="unknown", action = "write")
        write ( 205, '(A)' ) "# 1=#V, 2=|r|/|r0|, 3=|r_u|/|r_u0|, 4=|r_p|/|r_p0|"
        write ( 205, '(A)' ) "# 5=|Pr|/|Pr0|, 6=|Pr_u|/|Pr_u0|, 7=|Pr_p|/|Pr_p0|, 8=|Pr_est|/|Pr|"
      end if
    end if 

    nlevs=mla%nlevel
    dm= mla%dim

    do n=1,nlevs
       do d=1,dm
          call multifab_build_edge(r_u(n,d)  ,mla%la(n),1                ,1,d)
          call multifab_build_edge(w_u(n,d)  ,mla%la(n),1                ,0,d)
          call multifab_build_edge(tmp_u(n,d),mla%la(n),1                ,0,d)
          call multifab_build_edge(V_u(n,d)  ,mla%la(n),gmres_max_inner+1,0,d) ! Krylov vectors
       end do
       call multifab_build(r_p(n)  ,mla%la(n),1                ,1)
       call multifab_build(w_p(n)  ,mla%la(n),1                ,0)
       call multifab_build(tmp_p(n),mla%la(n),1                ,0)
       call multifab_build(V_p(n)  ,mla%la(n),gmres_max_inner+1,0) ! Krylov vectors
    end do

    ! apply scaling factor
    if (scale_factor .ne. 1.d0) then
       theta_alpha = theta_alpha*scale_factor
       do n=1,nlevs
          ! we will solve for scale*x_p so we need to scale the initial guess
          call multifab_mult_mult_s(x_p(n),scale_factor,x_p(n)%ng)
          ! scale the rhs:
          do i=1,dm
             call multifab_mult_mult_s(b_u(n,i),scale_factor,b_u(n,i)%ng)
          end do          
          ! scale the viscosities:
          call multifab_mult_mult_s(beta(n),scale_factor,beta(n)%ng)
          call multifab_mult_mult_s(gamma(n),scale_factor,gamma(n)%ng)
          do i=1,size(beta_ed,dim=2)
             call multifab_mult_mult_s(beta_ed(n,i),scale_factor,beta_ed(n,i)%ng)
          end do
       end do
    end if

    ! preconditioned norm_b: norm_pre_b
    call apply_precon(mla,b_u,b_p,tmp_u,tmp_p,alpha_fc, &
                      beta,beta_ed,gamma,theta_alpha,dx,the_bc_tower)
    call stag_l2_norm(mla,tmp_u,norm_u)
    call cc_l2_norm(mla,tmp_p,norm_p)
    norm_p=p_norm_weight*norm_p
    norm_pre_b = sqrt(norm_u**2+norm_p**2)

    if (present(norm_pre_rhs)) then
       norm_pre_rhs = norm_pre_b
    end if

    ! calculate the l2 norm of rhs
    call stag_l2_norm(mla,b_u,norm_u)
    call cc_l2_norm(mla,b_p,norm_p)
    norm_p = p_norm_weight*norm_p
    norm_b = sqrt(norm_u**2+norm_p**2)
    
    ! If norm_b=0 we should return zero as the solution and "return" from this routine
    if (norm_b .le. 0.0d0) then
       do n=1,nlevs
         do d=1,dm
           call setval(x_u(n,d),0.d0,all=.true.)
         end do 
         call setval(x_p(n),0.d0,all=.true.)
       end do   
       if (gmres_verbose .ge. 1) then
          if (parallel_IOProcessor()) then 
             write ( *, *) "gmres.f90: conveged in 0 iteration since rhs=0"
          end if 
       end if

       ! clean up memory
       do n=1,nlevs
          do d=1,dm
             call multifab_destroy(r_u(n,d))
             call multifab_destroy(w_u(n,d))
             call multifab_destroy(tmp_u(n,d))
             call multifab_destroy(V_u(n,d))
          end do
          call multifab_destroy(r_p(n))
          call multifab_destroy(w_p(n))
          call multifab_destroy(tmp_p(n))
          call multifab_destroy(V_p(n))
       end do

       return
    end if

    !!!!!!!!!!!!!!!!!!!
    ! begin outer iteration
    !!!!!!!!!!!!!!!!!!!
    
    total_iter = 0
    iter = 0  
    vcycle_counter = 0
    OuterLoop: do

       ! Calculate tmp = Ax
       call apply_matrix(mla,tmp_u,tmp_p,x_u,x_p,alpha_fc,beta,beta_ed, &
                         gamma,theta_alpha,dx,the_bc_tower)

       ! tmp = b - Ax
       do n=1,nlevs
          do d=1,dm
             call multifab_sub_sub_c(tmp_u(n,d),1,b_u(n,d),1,1,0)
             call multifab_mult_mult_s_c(tmp_u(n,d),1,-1.d0,1,0)
          end do
          call multifab_sub_sub_c(tmp_p(n),1,b_p(n),1,1,0)
          call multifab_mult_mult_s_c(tmp_p(n),1,-1.d0,1,0)
       end do
       
       ! un-preconditioned residuals
       call stag_l2_norm(mla,tmp_u, norm_u_noprecon)
       call cc_l2_norm(mla,tmp_p,norm_p_noprecon)
       norm_p_noprecon = p_norm_weight*norm_p_noprecon
       norm_resid_Stokes=sqrt(norm_u_noprecon**2+norm_p_noprecon**2)
       if(iter==0) then
         norm_init_u=norm_u_noprecon
         norm_init_p=norm_p_noprecon
         norm_init_Stokes=norm_resid_Stokes               
       end if 

       if (gmres_verbose .ge. 2) then
          if (parallel_IOProcessor()) then
             write ( *, '(a,i4)' )     'total ITERs = ', total_iter 
             write ( *, '(a,100g17.9)' )  'r/(r_0,b) = ',  norm_resid_Stokes/norm_init_Stokes, norm_resid_Stokes/norm_b
          end if
       end if
       if (gmres_verbose .ge. 3) then  
          if (parallel_IOProcessor()) then
             write ( *, '(a,100g17.9)' )  'un-Precond. rel. resid. (u,v,p) = ', &
                  norm_resid_Stokes/norm_init_Stokes, norm_u_noprecon/norm_init_Stokes, norm_p_noprecon/norm_init_Stokes
          end if
       end if
    
       ! solve for r = M^{-1} tmp
       ! We should not be counting these toward the number of mg cycles performed
       vcycle_counter_temp = vcycle_counter
       call apply_precon(mla,tmp_u,tmp_p,r_u,r_p,alpha_fc, &
                         beta,beta_ed,gamma,theta_alpha,dx,the_bc_tower)
       vcycle_counter = vcycle_counter_temp

       ! resid = sqrt(dot_product(r, r))
       call stag_l2_norm(mla,r_u,norm_u)
       call cc_l2_norm(mla,r_p,norm_p)
       norm_p = p_norm_weight*norm_p
       norm_resid=sqrt(norm_u**2+norm_p**2)
       ! If first iteration, save the initial preconditioned residual
       if(iter==0) then
         norm_init_resid=norm_resid
         norm_init_resid_u=norm_u
         norm_init_resid_p=norm_p
         norm_resid_est=norm_resid
       end if  

       if (gmres_verbose .ge. 3) then
         if (parallel_IOProcessor()) then 
           ! Write the residuals in outer iteration to 205
           write ( 205, '(I0,100g17.9)' ) vcycle_counter, &
              norm_resid/norm_init_resid, norm_u/norm_init_resid, norm_p/norm_init_resid, &
              norm_resid_Stokes/norm_init_Stokes, norm_u_noprecon/norm_init_Stokes, norm_p_noprecon/norm_init_Stokes, &
              norm_resid_est/norm_resid
           write ( *, '(a,100g17.9)' ) 'Precond. rel. res. (u,v,p) = ', &
                 norm_resid/norm_init_resid, norm_u/norm_init_resid, norm_p/norm_init_resid
                 ! Printing norm_u/norm_init_resid_u, norm_p/norm_init_resid_p may divide by zero...
         end if 
       end if
    
       ! We need to test the residual now and exit OuterLoop if converged
       if(total_iter >= gmres_max_iter) then
         if (gmres_verbose .ge. 1) then
           if (parallel_IOProcessor()) then
             write ( *, '(a)' ) 'GMRES did not converge in max number of total inner iterations: Exiting'
           end if
         end if
         exit OuterLoop
       else if(total_iter .ge. gmres_min_iter) then 
          ! other options
          if(norm_resid <= gmres_rel_tol*min(norm_pre_b, norm_init_resid)) then
          !if(norm_resid <= gmres_rel_tol*norm_pre_b) then
          !if(norm_resid <= gmres_rel_tol*norm_init_resid) then
            if (gmres_verbose .ge. 2) then
              if (parallel_IOProcessor()) then
                write ( *, '(a, i4,a,i4,a,i4)' ) 'GMRES converged: Outer = ', iter, ',  Inner = ', i, ' Total=', total_iter 
              end if
            end if  

            if (norm_resid_Stokes >= 10*gmres_rel_tol*min(norm_b, norm_init_Stokes)) then
              if (parallel_IOProcessor()) write(*,*) 'gmres.f90: Warning: gmres may not have converged: |r|/|b|=', &
                         norm_resid_Stokes/norm_b, ' |r|/|r0|=', norm_resid_Stokes/norm_init_Stokes
            end if

            exit OuterLoop ! Only exit if the *true* preconditioned residual is less than tolerance: Do not trust the gmres estimate
         else if (norm_resid <= gmres_abs_tol) then

            if (gmres_verbose .ge. 2) then
              if (parallel_IOProcessor()) then
                write ( *, '(a, i4,a,i4,a,i4)' ) 'GMRES converged: Outer = ', iter, ',  Inner = ', i, ' Total=', total_iter 
              end if
            end if  

            exit OuterLoop

          end if
       end if
       
       if(iter >= gmres_max_outer) then
           if (parallel_IOProcessor()) then 
               print*, 'GMRES did not converge in max number of outer iterations: Exiting'
           end if 
           exit OuterLoop       
       end if    
       iter=iter+1

       ! create the first basis in Krylov space
       ! V(1) = r / norm(r)
       do n=1,nlevs
          do d=1,dm
             call multifab_copy_c(V_u(n,d),1,r_u(n,d),1,1,0)
             call multifab_div_div_s_c(V_u(n,d),1,norm_resid,1,0)
          end do
          call multifab_copy_c(V_p(n),1,r_p(n),1,1,0)
          call multifab_div_div_s_c(V_p(n),1,norm_resid,1,0)
       end do

       ! s = norm(r) * e_1
       s(:) = 0.d0
       s(1) = norm_resid

       !!!!!!!!!!!!!!!!!!!!!!!!!
       ! begin inner iteration
       !!!!!!!!!!!!!!!!!!!!!!!!!

       InnerLoop: do i=1,gmres_max_inner
          total_iter = total_iter + 1
          i_copy=i

          ! tmp=A*V(i)
          ! we use r_p and r_u as temporaries to hold ith component of V
          do n=1,nlevs
             do d=1,dm
                call multifab_copy_c(r_u(n,d),1,V_u(n,d),i,1,0)
             end do
             call multifab_copy_c(r_p(n),1,V_p(n),i,1,0)
          end do
          
          call apply_matrix(mla,tmp_u,tmp_p,r_u,r_p,alpha_fc,beta,beta_ed, &
                            gamma,theta_alpha,dx,the_bc_tower)

          ! w = M^{-1} A*V(i)
          call apply_precon(mla,tmp_u,tmp_p,w_u,w_p,alpha_fc, &
                            beta,beta_ed,gamma,theta_alpha,dx,the_bc_tower)

          do k=1,i
             ! form H(k,i) Hessenberg matrix
             ! H(k,i) = dot_product(w, V(k))
             !        = dot_product(w_u, V_u(k))+dot_product(w_p, V_p(k))
             call stag_inner_prod(mla,w_u,1,V_u,k,inner_prod_vel)
             call cc_inner_prod(mla,w_p,1,V_p,k,inner_prod_pres)
             H(k,i) = sum(inner_prod_vel(1:dm)) + (p_norm_weight**2)*inner_prod_pres

             ! w = w - H(k,i) * V(k)
             ! use tmp_u and tmp_p as temporaries to hold kth component of V(k)
             do n=1,nlevs
                do d=1,dm
                   call multifab_copy_c(tmp_u(n,d),1,V_u(n,d),k,1,0)
                   call multifab_mult_mult_s_c(tmp_u(n,d),1,H(k,i),1,0)
                   call multifab_sub_sub_c(w_u(n,d),1,tmp_u(n,d),1,1,0)
                end do
                call multifab_copy_c(tmp_p(n),1,V_p(n),k,1,0)
                call multifab_mult_mult_s_c(tmp_p(n),1,H(k,i),1,0)
                call multifab_sub_sub_c(w_p(n),1,tmp_p(n),1,1,0)
             end do
          end do

          ! H(i+1,i) = norm(w)
          call stag_l2_norm(mla,w_u,norm_u)
          call cc_l2_norm(mla,w_p,norm_p)
          norm_p = p_norm_weight*norm_p
          H(i+1,i) = sqrt(norm_u**2+norm_p**2)

          ! V(i+1) = w / H(i+1,i)
          if ( H(i+1,i) .ne. 0.d0 ) then
             do n=1,nlevs
                do d=1,dm
                   call multifab_copy_c(V_u(n,d),i+1,w_u(n,d),1,1,0)
                   call multifab_div_div_s_c(V_u(n,d),i+1,H(i+1,i),1,0)
                end do
                call multifab_copy_c(V_p(n),i+1,w_p(n),1,1,0)
                call multifab_div_div_s_c(V_p(n),i+1,H(i+1,i),1,0)
             end do
          else 
             call bl_error("gmres.f90: error in orthogonalization")
          end if

          call least_squares(i, H, cs, sn, s)   ! solve least square problem
          norm_resid_est=abs(s(i+1)) 

          if (gmres_verbose .ge. 2) then
          if (parallel_IOProcessor()) then 
             write ( *, '(2i6,a,100g17.9)' ) total_iter, vcycle_counter, ',  est. rel. resid. |Pr|/(Pr0,b)= ', &
                norm_resid_est/norm_init_resid, norm_resid_est/norm_pre_b
             if (gmres_verbose .ge. 3) then
                ! Write the convergence information to files
                write ( 204, '(I0,100g17.9)' ) vcycle_counter, norm_resid_est/norm_init_resid
             end if
          end if 
          end if

          ! force output to be written for screen
          ! useful for parallel problems
          ! "6" is the standard output
          call flush(6)

          if(total_iter >= gmres_max_iter) then
            exit InnerLoop
          else if(total_iter .ge. gmres_min_iter) then 
             ! other options
             if ((norm_resid_est <= gmres_rel_tol*min(norm_pre_b, norm_init_resid)) .or. (norm_resid_est <= gmres_abs_tol)) then
             !if(norm_resid_est <= gmres_rel_tol*norm_pre_b .or. (norm_resid_est <= gmres_abs_tol)) then
             !if(norm_resid_est <= gmres_rel_tol*norm_init_resid .or. (norm_resid_est <= gmres_abs_tol)) then
                exit InnerLoop
             end if
          end if

       end do InnerLoop ! end of inner loop (do i=1,gmres_max_inner)

       ! update the solution   ! first, solve for y
       call SolveUTriangular(i_copy-1, H, s, y)

       ! then, x = x + dot(V(1:i),y(1:i))
       call update_sol(mla, x_u(1,:), x_p(1), V_u(1,:), V_p(1), y, i_copy, &
                       the_bc_tower%bc_tower_array(1), dx(1,:))

    end do OuterLoop ! end of outer loop (do iter=1,gmres_max_outer)
 
    ! AJN - this is here since I notice epsilon roundoff errors building up
    !       just enough to destroy the asymmetry in time-advancement codes that
    !       ultimately causes lack of convergence in subsequent gmres calls
    do i=1,dm
       call multifab_internal_sync(x_u(1,i))
    end do

    ! apply scaling factor
    if (scale_factor .ne. 1.d0) then
       theta_alpha = theta_alpha/scale_factor
       do n=1,nlevs
          ! the solution we got is scale*x_p
          call multifab_mult_mult_s(x_p(n),1.d0/scale_factor,x_p(n)%ng)
          ! unscale the rhs
          do i=1,dm
             call multifab_mult_mult_s(b_u(n,i),1.d0/scale_factor,b_u(n,i)%ng)
          end do
          ! unscale the viscosities
          call multifab_mult_mult_s(beta(n),1.d0/scale_factor,beta(n)%ng)
          call multifab_mult_mult_s(gamma(n),1.d0/scale_factor,gamma(n)%ng)
          do i=1,size(beta_ed,dim=2)
             call multifab_mult_mult_s(beta_ed(n,i),1.d0/scale_factor,beta_ed(n,i)%ng)
          end do
       end do
    end if

    ! clean up memory
    do n=1,nlevs
       do d=1,dm
          call multifab_destroy(r_u(n,d))
          call multifab_destroy(w_u(n,d))
          call multifab_destroy(tmp_u(n,d))
          call multifab_destroy(V_u(n,d))
       end do
       call multifab_destroy(r_p(n))
       call multifab_destroy(w_p(n))
       call multifab_destroy(tmp_p(n))
       call multifab_destroy(V_p(n))
    end do

    if (gmres_verbose .ge. 1) then
      if (gmres_verbose .ge. 3) then
        if (parallel_IOProcessor()) then 
          close(204)
          close(205)
        end if
      else
        if (parallel_IOProcessor()) then 
          write ( *, '(a)' ) 'Preconditioned GMRES:'
          write ( *, '(a,i4)' ) '  total ITERs = ', total_iter 
          write ( *, '(a, 100g17.9)' ) '  residual/(norm_b,initial) = ', norm_resid/norm_b, norm_resid/norm_init_resid 
        end if 
      end if
    end if

    call destroy(bpt)

  contains 

    !------------------------------------------
    subroutine update_sol(mla,x_u,x_p,V_u,V_p,y,i,the_bc_level,dx)

      type(ml_layout), intent(in   ) :: mla
      type(multifab) , intent(inout) :: x_u(:)
      type(multifab) , intent(inout) :: x_p
      type(multifab) , intent(inout) :: V_u(:)
      type(multifab) , intent(inout) :: V_p
      real(kind=dp_t), intent(in   ) :: y(:)
      integer        , intent(in   ) :: i
      type(bc_level) , intent(in   ) :: the_bc_level
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local variables
      integer :: dm,d,iter

      type(bl_prof_timer), save :: bpt

      call build(bpt,"udate_sol")

      dm = mla%dim

      ! set V(i) = V(i)*y(i)
      ! set x = x + V(i)
      do iter=1,i
         call multifab_mult_mult_s_c(V_p,iter,y(iter),1,0)
         call multifab_plus_plus_c(x_p,1,V_p,iter,1,0)
         do d=1,dm
            call multifab_mult_mult_s_c(V_u(d),iter,y(iter),1,0)
            call multifab_plus_plus_c(x_u(d),1,V_u(d),iter,1,0)
         end do
      end do
      
      call destroy(bpt)

    end subroutine update_sol

    subroutine least_squares(i, H, cs, sn, s)

      integer        , intent(in   ) :: i
      real(kind=dp_t), intent(inout) :: H(:,:)
      real(kind=dp_t), intent(inout) :: cs(:)
      real(kind=dp_t), intent(inout) :: sn(:)
      real(kind=dp_t), intent(inout) :: s(:)

      ! local variable       
      integer :: k
      real(kind=dp_t) :: temp

      type(bl_prof_timer), save :: bpt      

      call build(bpt,"least_squares")

      ! apply Givens rotation
      do k = 1,i-1
         temp     =  cs(k)*H(k,i) + sn(k)*H(k+1,i)
         H(k+1,i) = -sn(k)*H(k,i) + cs(k)*H(k+1,i)
         H(k,i)   = temp
      end do

      ! form i-th rotation matrix
      call rotmat(H(i,i), H(i+1,i), cs(i), sn(i))

      ! approximate residual norm
      temp   = cs(i)*s(i)                        
      s(i+1) = -sn(i)*s(i)
      s(i)   = temp
      H(i,i) = cs(i)*H(i,i) + sn(i)*H(i+1,i)
      H(i+1,i) = 0.0d0

      call destroy(bpt)

    end subroutine least_squares

    subroutine rotmat(a, b, cs, sn)

      ! Compute the Givens rotation matrix parameters for a and b.

      real(kind=dp_t), intent(in   ) :: a,b
      real(kind=dp_t), intent(inout) :: cs,sn

      ! local
      real(kind=dp_t) :: temp

      type(bl_prof_timer), save :: bpt

      call build(bpt,"rotmat")

      if ( b .eq. 0.0d0 ) then
         cs = 1.0d0
         sn = 0.0d0
      elseif(abs(b) > abs(a)) then
         temp = a/b
         sn = 1.0d0/sqrt(1.0d0 + temp**2)
         cs = temp * sn
      else
         temp = b / a;
         cs = 1.0d0/sqrt(1.0d0 + temp**2)
         sn = temp * cs
      end if

      call destroy(bpt)

    end subroutine rotmat

    subroutine SolveUTriangular(k, H, s, y)

      integer        , intent(in   ) :: k
      real(kind=dp_t), intent(in   ) :: H(:,:)
      real(kind=dp_t), intent(in   ) :: s(:)
      real(kind=dp_t), intent(inout) :: y(:)
      !local 
      integer :: i

      type(bl_prof_timer), save :: bpt

      call build(bpt,"SolveUTriangular")

      y(k+1) = s(k+1)/H(k+1,k+1)
      do i = k, 1, -1
         y(i) = (s(i)-dot_product(H(i,i+1:k+1), y(i+1:k+1)))/H(i,i)
      end do

      call destroy(bpt)

    end subroutine SolveUTriangular

  end subroutine gmres

end module gmres_module
