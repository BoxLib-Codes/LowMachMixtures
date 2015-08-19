module stag_applyop_module

  use ml_layout_module
  use multifab_module
  use probin_common_module, only: visc_type
  use define_bc_module
  use bc_module
  use multifab_physbc_stag_module
  use convert_stag_module

  implicit none

  private

  public :: stag_applyop, stag_applyop_level

contains
  
  ! compute Lphi
  subroutine stag_applyop(mla,the_bc_tower,phi_fc,Lphi_fc,alpha_fc, &
                          beta_cc,beta_ed,gamma_cc,theta_alpha,dx)

    type(ml_layout), intent(in   ) :: mla
    type(bc_tower) , intent(in   ) :: the_bc_tower
    type(multifab) , intent(in   ) :: phi_fc(:,:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:,:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_fc(:,:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc(:)    ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:,:)  ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma_cc(:)   ! cell-centered
    real(kind=dp_t), intent(in   ) :: theta_alpha,dx(:,:)

    ! local
    integer :: i,n,dm,nlevs
    type(multifab) :: alpha_fc_temp(mla%nlevel,mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"stag_applyop")

    dm = mla%dim
    nlevs = mla%nlevel

    ! multiply alpha_fc_temp by theta_alpha
    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(alpha_fc_temp(n,i),mla%la(n),1,0,i)
          call multifab_copy_c(alpha_fc_temp(n,i),1,alpha_fc(n,i),1,1,0)
          call multifab_mult_mult_s_c(alpha_fc_temp(n,i),1,theta_alpha,1,0)
       end do
    end do

    do n=1,nlevs
       call stag_applyop_level(mla%la(n),the_bc_tower%bc_tower_array(n), &
                               phi_fc(n,:),Lphi_fc(n,:),alpha_fc_temp(n,:), &
                               beta_cc(n),beta_ed(n,:),gamma_cc(n),dx(n,:))
    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(alpha_fc_temp(n,i))
       end do
    end do
       
    call destroy(bpt)

  end subroutine stag_applyop

  ! compute Lphi
  subroutine stag_applyop_level(la,the_bc_level,phi_fc,Lphi_fc,alpha_fc,beta_cc, &
                                beta_ed,gamma_cc,dx,color_in)
    
    type(layout)   , intent(in   ) :: la
    type(bc_level) , intent(in   ) :: the_bc_level
    type(multifab) , intent(in   ) :: phi_fc(:)   ! face-centered
    type(multifab) , intent(inout) :: Lphi_fc(:)  ! face-centered
    type(multifab) , intent(in   ) :: alpha_fc(:) ! face-centered
    type(multifab) , intent(in   ) :: beta_cc     ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:)  ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(in   ) :: gamma_cc    ! cell-centered
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ), optional :: color_in
    
    ! local
    integer :: i,dm,ng_p,ng_l,ng_a,ng_b,ng_e,ng_g
    integer :: lo(get_dim(la)), hi(get_dim(la))
    integer :: color
    
    real(kind=dp_t), pointer :: ppx(:,:,:,:)
    real(kind=dp_t), pointer :: ppy(:,:,:,:)
    real(kind=dp_t), pointer :: ppz(:,:,:,:)
    real(kind=dp_t), pointer :: lpx(:,:,:,:)
    real(kind=dp_t), pointer :: lpy(:,:,:,:)
    real(kind=dp_t), pointer :: lpz(:,:,:,:)
    real(kind=dp_t), pointer :: apx(:,:,:,:)
    real(kind=dp_t), pointer :: apy(:,:,:,:)
    real(kind=dp_t), pointer :: apz(:,:,:,:)
    real(kind=dp_t), pointer ::  bp(:,:,:,:)
    real(kind=dp_t), pointer :: bp1(:,:,:,:)
    real(kind=dp_t), pointer :: bp2(:,:,:,:)
    real(kind=dp_t), pointer :: bp3(:,:,:,:)
    real(kind=dp_t), pointer :: bnp(:,:,:,:)
    real(kind=dp_t), pointer :: kp(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(beta_cc%dim), xhi(beta_cc%dim)
    integer :: ylo(beta_cc%dim), yhi(beta_cc%dim)
    integer :: zlo(beta_cc%dim), zhi(beta_cc%dim)


    type(bl_prof_timer),save :: bpt

    call build(bpt,"stag_applyop_level")

    if (dx(1) .ne. dx(2)) then
       call bl_error("stag_applyop_2d requires the same dx in all directions")
    end if

    dm = get_dim(la)

    if (present(color_in)) then
       color = color_in
    else
       color = 0
    end if

    ng_p = phi_fc(1)%ng
    ng_l = Lphi_fc(1)%ng
    ng_a = alpha_fc(1)%ng
    ng_b = beta_cc%ng
    ng_e = beta_ed(1)%ng
    ng_g = gamma_cc%ng

    !$omp parallel private(mfi,i,xnodalbox,ynodalbox,znodalbox,xlo,ylo,zlo) &
    !$omp private(xhi,yhi,zhi,ppx,ppy,ppz,lpx,lpy,lpz,apx,apy,apz,bp,kp,lo,hi) &
    !$omp private(bp1,bp2,bp3)

    call mfiter_build(mfi, beta_cc, tiling=.true.)

    do while (more_tile(mfi))
       i = get_fab_index(mfi)

       xnodalbox = get_nodaltilebox(mfi,1)
       xlo = lwb(xnodalbox)
       xhi = upb(xnodalbox)
       ynodalbox = get_nodaltilebox(mfi,2)
       ylo = lwb(ynodalbox)
       yhi = upb(ynodalbox)
       znodalbox = get_nodaltilebox(mfi,3)
       zlo = lwb(znodalbox)
       zhi = upb(znodalbox)

!    do i=1,nfabs(phi_fc(1))
       ppx => dataptr(phi_fc(1), i)
       ppy => dataptr(phi_fc(2), i)
       lpx => dataptr(Lphi_fc(1), i)
       lpy => dataptr(Lphi_fc(2), i)
       apx => dataptr(alpha_fc(1), i)
       apy => dataptr(alpha_fc(2), i)
       bp  => dataptr(beta_cc, i)
       kp  => dataptr(gamma_cc, i)
       lo = lwb(get_box(phi_fc(1), i))
       hi = upb(get_box(phi_fc(1), i))
       select case (dm)
       case (2)
          bnp => dataptr(beta_ed(1), i)
          call stag_applyop_2d(ppx(:,:,1,1),ppy(:,:,1,1),ng_p, &
                               lpx(:,:,1,1),lpy(:,:,1,1),ng_l, &
                               apx(:,:,1,1),apy(:,:,1,1),ng_a, &
                               bp(:,:,1,1),ng_b, &
                               bnp(:,:,1,1),ng_e, &
                               kp(:,:,1,1),ng_g, &
                               lo,hi,dx,color,xlo,xhi,ylo,yhi)
       case (3)
       ppz => dataptr(phi_fc(3), i)
       lpz => dataptr(Lphi_fc(3), i)
       apz => dataptr(alpha_fc(3), i)
       bp1 => dataptr(beta_ed(1), i)
       bp2 => dataptr(beta_ed(2), i)
       bp3 => dataptr(beta_ed(3), i)
       call stag_applyop_3d(ppx(:,:,:,1),ppy(:,:,:,1),ppz(:,:,:,1),ng_p, &
                            lpx(:,:,:,1),lpy(:,:,:,1),lpz(:,:,:,1),ng_l, &
                            apx(:,:,:,1),apy(:,:,:,1),apz(:,:,:,1),ng_a, &
                            bp(:,:,:,1),ng_b, &
                            bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                            kp(:,:,:,1),ng_g, &
                            lo,hi,dx,color,xlo,xhi,ylo,yhi,zlo,zhi)

       end select
    end do

    !$omp end parallel

    do i=1,dm
       ! set Lphi on physical domain boundaries to zero
       call multifab_physbc_domainvel(Lphi_fc(i),vel_bc_comp+i-1,the_bc_level,dx)
    end do

    call destroy(bpt)

  end subroutine stag_applyop_level

  subroutine stag_applyop_2d(phix,phiy,ng_p,Lpx,Lpy,ng_l, &
                             alphax,alphay,ng_a,beta,ng_b,beta_ed,ng_n, &
                             gamma,ng_g,glo,ghi,dx,color,xlo,xhi,ylo,yhi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_p,ng_l,ng_a,ng_b,ng_n,ng_g
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
    real(kind=dp_t), intent(in   ) ::    phix(glo(1)-ng_p:,glo(2)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiy(glo(1)-ng_p:,glo(2)-ng_p:)
    real(kind=dp_t), intent(inout) ::     Lpx(glo(1)-ng_l:,glo(2)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpy(glo(1)-ng_l:,glo(2)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(glo(1)-ng_a:,glo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(glo(1)-ng_a:,glo(2)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(glo(1)-ng_b:,glo(2)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_ed(glo(1)-ng_n:,glo(2)-ng_n:)
    real(kind=dp_t), intent(in   ) ::   gamma(glo(1)-ng_g:,glo(2)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: color
    
    ! local
    integer :: i,j

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds, dxsqinv
    real(kind=dp_t) :: b,c

    ! coloring parameters
    logical :: do_x, do_y
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       offset = 2
    end if

    dxsq = dx(1)**2
    dxsqinv = 1.d0/dxsq
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0

    if (visc_type .eq. -1) then

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     (beta(i,j)+beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))*dxsqinv) &
                     - ( phix(i+1,j)*beta(i,j) &
                     +phix(i-1,j)*beta(i-1,j) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) )*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     (beta(i,j)+beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))*dxsqinv) &
                     - ( phiy(i,j+1)*beta(i,j) &
                     +phiy(i,j-1)*beta(i,j-1) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) )*dxsqinv

             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(xlo(1),xlo(2))

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + 4.d0*b*dxsqinv) &
                     -(phix(i+1,j)+phix(i-1,j)+phix(i,j+1)+phix(i,j-1))*b*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 4.d0*b*dxsqinv) &
                     -(phiy(i,j+1)+phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j))*b*dxsqinv

             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))*dxsqinv) &
                     
                     -( 2.d0*phix(i+1,j)*beta(i,j) &
                     +2.d0*phix(i-1,j)*beta(i-1,j) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) &
                     
                     +phiy(i,j+1)*beta_ed(i,j+1) &
                     -phiy(i,j)*beta_ed(i,j) &
                     -phiy(i-1,j+1)*beta_ed(i,j+1) &
                     +phiy(i-1,j)*beta_ed(i,j) )*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))*dxsqinv) &
                     
                     -( 2.d0*phiy(i,j+1)*beta(i,j) &
                     +2.d0*phiy(i,j-1)*beta(i,j-1) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) &
                     
                     +phix(i+1,j)*beta_ed(i+1,j) &
                     -phix(i,j)*beta_ed(i,j) &
                     -phix(i+1,j-1)*beta_ed(i+1,j) &
                     +phix(i,j-1)*beta_ed(i,j) )*dxsqinv

             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(xlo(1),xlo(2))

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + 6.d0*b*dxsqinv) &
                     -(2.d0*phix(i+1,j)+2.d0*phix(i-1,j)+phix(i,j+1)+phix(i,j-1) &
                     +phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*b*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + 6.d0*b*dxsqinv) &
                     -(2.d0*phiy(i,j+1)+2.d0*phiy(i,j-1)+phiy(i+1,j)+phiy(i-1,j) &
                     +phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*b*dxsqinv

             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i-1,j)+gamma(i-1,j) &
                     +beta_ed(i,j+1)+beta_ed(i,j) )*dxsqinv) &
                     
                     -( phix(i+1,j)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phix(i-1,j)*(fourthirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phix(i,j+1)*beta_ed(i,j+1) &
                     +phix(i,j-1)*beta_ed(i,j) &
                     
                     +phiy(i,j+1)*(beta_ed(i,j+1)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i,j)*(beta_ed(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phiy(i-1,j+1)*(beta_ed(i,j+1)-twothirds*beta(i-1,j)+gamma(i-1,j)) &
                     +phiy(i-1,j)*(beta_ed(i,j)-twothirds*beta(i-1,j)+gamma(i-1,j)) )*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j) + &
                     ( fourthirds*beta(i,j)+gamma(i,j)+fourthirds*beta(i,j-1)+gamma(i,j-1) &
                     +beta_ed(i+1,j)+beta_ed(i,j) )*dxsqinv) &
                     
                     -( phiy(i,j+1)*(fourthirds*beta(i,j)+gamma(i,j)) &
                     +phiy(i,j-1)*(fourthirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phiy(i+1,j)*beta_ed(i+1,j) &
                     +phiy(i-1,j)*beta_ed(i,j) &
                     
                     +phix(i+1,j)*(beta_ed(i+1,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i,j)*(beta_ed(i,j)-twothirds*beta(i,j)+gamma(i,j)) &
                     -phix(i+1,j-1)*(beta_ed(i+1,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) &
                     +phix(i,j-1)*(beta_ed(i,j)-twothirds*beta(i,j-1)+gamma(i,j-1)) )*dxsqinv

             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(xlo(1),xlo(2))
       c = gamma(xlo(1),xlo(2))

       if (do_x) then

          do j=xlo(2),xhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(xlo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=xlo(1)+ioff,xhi(1),offset

                Lpx(i,j) = phix(i,j)*(alphax(i,j)+(14.d0*b*onethird+2.d0*c)*dxsqinv) &
                     -((phix(i+1,j)+phix(i-1,j))*(fourthirds*b+c) &
                     +(phix(i,j+1)+phix(i,j-1))*b &
                     +(phiy(i,j+1)-phiy(i,j)-phiy(i-1,j+1)+phiy(i-1,j))*(onethird*b+c))*dxsqinv

             end do
          end do

       end if

       if (do_y) then

          do j=ylo(2),yhi(2)
             ioff = 0
             if ( offset .eq. 2 .and. mod(ylo(1)+j,2) .ne. mod(color+1,2) ) ioff = 1
             do i=ylo(1)+ioff,yhi(1),offset

                Lpy(i,j) = phiy(i,j)*(alphay(i,j)+(14.d0*b*onethird+2.d0*c)*dxsqinv) &
                     -((phiy(i,j+1)+phiy(i,j-1))*(fourthirds*b+c) &
                     +(phiy(i+1,j)+phiy(i-1,j))*b &
                     +(phix(i+1,j)-phix(i,j)-phix(i+1,j-1)+phix(i,j-1))*(onethird*b+c))*dxsqinv

             end do
          end do

       end if

    end if

  end subroutine stag_applyop_2d

  subroutine stag_applyop_3d(phix,phiy,phiz,ng_p,Lpx,Lpy,Lpz,ng_l, &
                             alphax,alphay,alphaz,ng_a,beta,ng_b, &
                             beta_xy,beta_xz,beta_yz,ng_e, &
                             gamma,ng_g,glo,ghi,dx,color,xlo,xhi,ylo,yhi,zlo,zhi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_p,ng_l,ng_a,ng_b,ng_e,ng_g
    integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t), intent(in   ) ::    phix(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiy(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:)
    real(kind=dp_t), intent(in   ) ::    phiz(glo(1)-ng_p:,glo(2)-ng_p:,glo(3)-ng_p:)
    real(kind=dp_t), intent(inout) ::     Lpx(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpy(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:)
    real(kind=dp_t), intent(inout) ::     Lpz(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:)
    real(kind=dp_t), intent(in   ) ::  alphax(glo(1)-ng_a:,glo(2)-ng_a:,glo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphay(glo(1)-ng_a:,glo(2)-ng_a:,glo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::  alphaz(glo(1)-ng_a:,glo(2)-ng_a:,glo(3)-ng_a:)
    real(kind=dp_t), intent(in   ) ::    beta(glo(1)-ng_b:,glo(2)-ng_b:,glo(3)-ng_b:)
    real(kind=dp_t), intent(in   ) :: beta_xy(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_xz(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_yz(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) ::   gamma(glo(1)-ng_g:,glo(2)-ng_g:,glo(3)-ng_g:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: color

    ! local
    integer :: i,j,k

    real(kind=dp_t) :: dxsq, onethird, twothirds, fourthirds,dxsqinv
    real(kind=dp_t) :: b,c

    ! coloring parameters
    logical :: do_x, do_y, do_z
    integer :: offset, ioff

    do_x = .true.
    do_y = .true.
    do_z = .true.
    offset = 1

    if (color .eq. 1 .or. color .eq. 2) then
       do_y = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 3 .or. color .eq. 4) then
       do_x = .false.
       do_z = .false.
       offset = 2
    else if (color .eq. 5 .or. color .eq. 6) then
       do_x = .false.
       do_y = .false.
       offset = 2
    end if

    dxsq = dx(1)**2
    dxsqinv = 1.d0/dxsq
    onethird = 1.d0/3.d0
    twothirds = 2.d0/3.d0
    fourthirds = 4.d0/3.d0
    
    if (visc_type .eq. -1) then

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + &
                        ( beta(i,j,k)+beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )*dxsqinv) &
                        - ( phix(i+1,j,k)*beta(i,j,k) &
                        +phix(i-1,j,k)*beta(i-1,j,k) &
                        +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                        +phix(i,j-1,k)*beta_xy(i,j,k) &
                        +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                        +phix(i,j,k-1)*beta_xz(i,j,k) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )*dxsqinv) &
                        - ( phiy(i,j+1,k)*beta(i,j,k) &
                        +phiy(i,j-1,k)*beta(i,j-1,k) &
                        +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                        +phiy(i-1,j,k)*beta_xy(i,j,k) &
                        +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                        +phiy(i,j,k-1)*beta_yz(i,j,k) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + &
                        ( beta(i,j,k)+beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )*dxsqinv) &
                        - ( phiz(i,j,k+1)*beta(i,j,k) &
                        +phiz(i,j,k-1)*beta(i,j,k-1) &
                        +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                        +phiz(i-1,j,k)*beta_xz(i,j,k) &
                        +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                        +phiz(i,j-1,k)*beta_yz(i,j,k) )*dxsqinv


                end do
             end do
          end do

       end if

    else if (visc_type .eq. 1) then

       b = beta(xlo(1),xlo(2),xlo(3))

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 6.d0*b*dxsqinv) &
                        -( phix(i+1,j,k)+phix(i-1,j,k) &
                        +phix(i,j+1,k)+phix(i,j-1,k) &
                        +phix(i,j,k+1)+phix(i,j,k-1))*b*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 6.d0*b*dxsqinv) &
                        -( phiy(i+1,j,k)+phiy(i-1,j,k) &
                        +phiy(i,j+1,k)+phiy(i,j-1,k) &
                        +phiy(i,j,k+1)+phiy(i,j,k-1))*b*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 6.d0*b*dxsqinv) &
                        -( phiz(i+1,j,k)+phiz(i-1,j,k) &
                        +phiz(i,j+1,k)+phiz(i,j-1,k) &
                        +phiz(i,j,k+1)+phiz(i,j,k-1))*b*dxsqinv

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -2) then

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*( alphax(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )*dxsqinv ) &
                        
                        -( 2.d0*phix(i+1,j,k)*beta(i,j,k) &
                        +2.d0*phix(i-1,j,k)*beta(i-1,j,k) &
                        +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                        +phix(i,j-1,k)*beta_xy(i,j,k) &
                        +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                        +phix(i,j,k-1)*beta_xz(i,j,k) &
                        
                        +phiy(i,j+1,k)*beta_xy(i,j+1,k) &
                        -phiy(i,j,k)*beta_xy(i,j,k) &
                        -phiy(i-1,j+1,k)*beta_xy(i,j+1,k) &
                        +phiy(i-1,j,k)*beta_xy(i,j,k) &
                        
                        +phiz(i,j,k+1)*beta_xz(i,j,k+1) &
                        -phiz(i,j,k)*beta_xz(i,j,k) &
                        -phiz(i-1,j,k+1)*beta_xz(i,j,k+1) &
                        +phiz(i-1,j,k)*beta_xz(i,j,k) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*( alphay(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )*dxsqinv ) &
                        
                        -( 2.d0*phiy(i,j+1,k)*beta(i,j,k) &
                        +2.d0*phiy(i,j-1,k)*beta(i,j-1,k) &
                        +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                        +phiy(i-1,j,k)*beta_xy(i,j,k) &
                        +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                        +phiy(i,j,k-1)*beta_yz(i,j,k) &
                        
                        +phix(i+1,j,k)*beta_xy(i+1,j,k) &
                        -phix(i,j,k)*beta_xy(i,j,k) &
                        -phix(i+1,j-1,k)*beta_xy(i+1,j,k) &
                        +phix(i,j-1,k)*beta_xy(i,j,k) &
                        
                        +phiz(i,j,k+1)*beta_yz(i,j,k+1) &
                        -phiz(i,j,k)*beta_yz(i,j,k) &
                        -phiz(i,j-1,k+1)*beta_yz(i,j,k+1) &
                        +phiz(i,j-1,k)*beta_yz(i,j,k) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*( alphaz(i,j,k) + &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )*dxsqinv ) &
                        
                        -( 2.d0*phiz(i,j,k+1)*beta(i,j,k) &
                        +2.d0*phiz(i,j,k-1)*beta(i,j,k-1) &
                        +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                        +phiz(i-1,j,k)*beta_xz(i,j,k) &
                        +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                        +phiz(i,j-1,k)*beta_yz(i,j,k) &
                        
                        +phix(i+1,j,k)*beta_xz(i+1,j,k) &
                        -phix(i,j,k)*beta_xz(i,j,k) &
                        -phix(i+1,j,k-1)*beta_xz(i+1,j,k) &
                        +phix(i,j,k-1)*beta_xz(i,j,k) &
                        
                        +phiy(i,j+1,k)*beta_yz(i,j+1,k) &
                        -phiy(i,j,k)*beta_yz(i,j,k) &
                        -phiy(i,j+1,k-1)*beta_yz(i,j+1,k) &
                        +phiy(i,j,k-1)*beta_yz(i,j,k) )*dxsqinv

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 2) then

       b = beta(xlo(1),xlo(2),xlo(3))

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k) + 8.d0*b*dxsqinv) &
                        -( 2.d0*phix(i+1,j,k)+2.d0*phix(i-1,j,k) &
                        +phix(i,j+1,k)+phix(i,j-1,k) &
                        +phix(i,j,k+1)+phix(i,j,k-1) &
                        +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k) )*b*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k) + 8.d0*b*dxsqinv) &
                        -( 2.d0*phiy(i,j+1,k)+2.d0*phiy(i,j-1,k) &
                        +phiy(i+1,j,k)+phiy(i-1,j,k) &
                        +phiy(i,j,k+1)+phiy(i,j,k-1) &
                        +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k) )*b*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k) + 8.d0*b*dxsqinv) &
                        -( 2.d0*phiz(i,j,k+1)+2.d0*phiz(i,j,k-1) &
                        +phiz(i+1,j,k)+phiz(i-1,j,k) &
                        +phiz(i,j+1,k)+phiz(i,j-1,k) &
                        +phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j,k-1)+phix(i,j,k-1) &
                        +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i,j+1,k-1)+phiy(i,j,k-1) )*b*dxsqinv

                end do
             end do
          end do

       end if

    else if (visc_type .eq. -3) then

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*( alphax(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i-1,j,k)+gamma(i-1,j,k) &
                        +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                        +beta_xz(i,j,k)+beta_xz(i,j,k+1) )*dxsqinv ) &
                        
                        -( phix(i+1,j,k)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                        +phix(i-1,j,k)*(fourthirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                        +phix(i,j+1,k)*beta_xy(i,j+1,k) &
                        +phix(i,j-1,k)*beta_xy(i,j,k) &
                        +phix(i,j,k+1)*beta_xz(i,j,k+1) &
                        +phix(i,j,k-1)*beta_xz(i,j,k) &
                        
                        +phiy(i,j+1,k)*(beta_xy(i,j+1,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiy(i,j,k)*(beta_xy(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiy(i-1,j+1,k)*(beta_xy(i,j+1,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                        +phiy(i-1,j,k)*(beta_xy(i,j,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                        
                        +phiz(i,j,k+1)*(beta_xz(i,j,k+1)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiz(i,j,k)*(beta_xz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiz(i-1,j,k+1)*(beta_xz(i,j,k+1)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) &
                        +phiz(i-1,j,k)*(beta_xz(i,j,k)-twothirds*beta(i-1,j,k)+gamma(i-1,j,k)) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*( alphay(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j-1,k)+gamma(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )*dxsqinv ) &
                        
                        -( phiy(i,j+1,k)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                        +phiy(i,j-1,k)*(fourthirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                        +phiy(i+1,j,k)*beta_xy(i+1,j,k) &
                        +phiy(i-1,j,k)*beta_xy(i,j,k) &
                        +phiy(i,j,k+1)*beta_yz(i,j,k+1) &
                        +phiy(i,j,k-1)*beta_yz(i,j,k) &
                        
                        +phix(i+1,j,k)*(beta_xy(i+1,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phix(i,j,k)*(beta_xy(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phix(i+1,j-1,k)*(beta_xy(i+1,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                        +phix(i,j-1,k)*(beta_xy(i,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                        
                        +phiz(i,j,k+1)*(beta_yz(i,j,k+1)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiz(i,j,k)*(beta_yz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiz(i,j-1,k+1)*(beta_yz(i,j,k+1)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) &
                        +phiz(i,j-1,k)*(beta_yz(i,j,k)-twothirds*beta(i,j-1,k)+gamma(i,j-1,k)) )*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*( alphaz(i,j,k) + &
                        ( fourthirds*beta(i,j,k)+gamma(i,j,k) &
                        +fourthirds*beta(i,j,k-1)+gamma(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )*dxsqinv ) &
                        
                        -( phiz(i,j,k+1)*(fourthirds*beta(i,j,k)+gamma(i,j,k)) &
                        +phiz(i,j,k-1)*(fourthirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                        +phiz(i+1,j,k)*beta_xz(i+1,j,k) &
                        +phiz(i-1,j,k)*beta_xz(i,j,k) &
                        +phiz(i,j+1,k)*beta_yz(i,j+1,k) &
                        +phiz(i,j-1,k)*beta_yz(i,j,k) &
                        
                        +phix(i+1,j,k)*(beta_xz(i+1,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phix(i,j,k)*(beta_xz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phix(i+1,j,k-1)*(beta_xz(i+1,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                        +phix(i,j,k-1)*(beta_xz(i,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                        
                        +phiy(i,j+1,k)*(beta_yz(i,j+1,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiy(i,j,k)*(beta_yz(i,j,k)-twothirds*beta(i,j,k)+gamma(i,j,k)) &
                        -phiy(i,j+1,k-1)*(beta_yz(i,j+1,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) &
                        +phiy(i,j,k-1)*(beta_yz(i,j,k)-twothirds*beta(i,j,k-1)+gamma(i,j,k-1)) )*dxsqinv

                end do
             end do
          end do

       end if

    else if (visc_type .eq. 3) then

       b = beta(xlo(1),xlo(2),xlo(3))
       c = gamma(xlo(1),xlo(2),xlo(3))

       if (do_x) then

          do k=xlo(3),xhi(3)
             do j=xlo(2),xhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(xlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=xlo(1)+ioff,xhi(1),offset

                   Lpx(i,j,k) = phix(i,j,k)*(alphax(i,j,k)+(20.d0*b*onethird+2.d0*c)*dxsqinv) &
                        - ((phix(i+1,j,k)+phix(i-1,j,k))*(fourthirds*b+c) &
                        +(phix(i,j+1,k)+phix(i,j-1,k)+phix(i,j,k+1)+phix(i,j,k-1))*b &
                        +( phiy(i,j+1,k)-phiy(i,j,k)-phiy(i-1,j+1,k)+phiy(i-1,j,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i-1,j,k+1)+phiz(i-1,j,k)) &
                        *(onethird*b+c))*dxsqinv

                end do
             end do
          end do

       end if

       if (do_y) then

          do k=ylo(3),yhi(3)
             do j=ylo(2),yhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(ylo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=ylo(1)+ioff,yhi(1),offset

                   Lpy(i,j,k) = phiy(i,j,k)*(alphay(i,j,k)+(20.d0*b*onethird+2.d0*c)*dxsqinv) &
                        - ((phiy(i,j+1,k)+phiy(i,j-1,k))*(fourthirds*b+c) &
                        +(phiy(i+1,j,k)+phiy(i-1,j,k)+phiy(i,j,k+1)+phiy(i,j,k-1))*b &
                        +( phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j-1,k)+phix(i,j-1,k) &
                        +phiz(i,j,k+1)-phiz(i,j,k)-phiz(i,j-1,k+1)+phiz(i,j-1,k)) &
                        *(onethird*b+c))*dxsqinv

                end do
             end do
          end do

       end if

       if (do_z) then

          do k=zlo(3),zhi(3)
             do j=zlo(2),zhi(2)
                ioff = 0
                if ( offset .eq. 2 .and. mod(zlo(1)+j+k,2) .ne. mod(color+1,2) ) ioff = 1
                do i=zlo(1)+ioff,zhi(1),offset

                   Lpz(i,j,k) = phiz(i,j,k)*(alphaz(i,j,k)+(20.d0*b*onethird+2.d0*c)*dxsqinv) &
                        - ((phiz(i,j,k+1)+phiz(i,j,k-1))*(fourthirds*b+c) &
                        +(phiz(i+1,j,k)+phiz(i-1,j,k)+phiz(i,j+1,k)+phiz(i,j-1,k))*b &
                        +( phix(i+1,j,k)-phix(i,j,k)-phix(i+1,j,k-1)+phix(i,j,k-1) &
                        +phiy(i,j+1,k)-phiy(i,j,k)-phiy(i,j+1,k-1)+phiy(i,j,k-1)) &
                        *(onethird*b+c))*dxsqinv

                end do
             end do
          end do

       end if

    end if

  end subroutine stag_applyop_3d

end module stag_applyop_module
