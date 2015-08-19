module inverse_diag_lap_module

  use ml_layout_module
  use multifab_module
  use probin_common_module, only: visc_type

  implicit none

  private

  public :: inverse_diag_lap

contains
  
  ! take a coefficient and compute the inverse diagonal of the weighted Laplacian operator
  ! ignoring the dx^2 scaling
  subroutine inverse_diag_lap(mla,beta_cc,beta_ed,Minv_fc)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: beta_cc(:)   ! cell-centered
    type(multifab) , intent(in   ) :: beta_ed(:,:) ! nodal (2d); edge-centered (3d)
    type(multifab) , intent(inout) :: Minv_fc(:,:) ! face-centered
    
    ! local
    integer :: n,nlevs,i,dm,ng_c,ng_e,ng_f
    integer :: lo(mla%dim),hi(mla%dim)
    
    real(kind=dp_t), pointer :: bpc(:,:,:,:)
    real(kind=dp_t), pointer :: bp1(:,:,:,:)
    real(kind=dp_t), pointer :: bp2(:,:,:,:)
    real(kind=dp_t), pointer :: bp3(:,:,:,:)
    real(kind=dp_t), pointer :: bnp(:,:,:,:)
    real(kind=dp_t), pointer :: Mpx(:,:,:,:)
    real(kind=dp_t), pointer :: Mpy(:,:,:,:)
    real(kind=dp_t), pointer :: Mpz(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"inverse_diag_lap")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_c = beta_cc(1)%ng
    ng_e = beta_ed(1,1)%ng
    ng_f = Minv_fc(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(beta_cc(n))
          Mpx => dataptr(Minv_fc(n,1), i)
          Mpy => dataptr(Minv_fc(n,2), i)
          bpc => dataptr(beta_cc(n), i)
          lo = lwb(get_box(beta_cc(n), i))
          hi = upb(get_box(beta_cc(n), i))
          select case (dm)
          case (2)
             bnp => dataptr(beta_ed(n,1), i)
             call inverse_diag_lap_2d(Mpx(:,:,1,1),Mpy(:,:,1,1),ng_f, &
                                      bpc(:,:,1,1),ng_c, &
                                      bnp(:,:,1,1),ng_e, lo, hi)
          case (3)
             bp1 => dataptr(beta_ed(n,1), i)
             bp2 => dataptr(beta_ed(n,2), i)
             bp3 => dataptr(beta_ed(n,3), i)
             call inverse_diag_lap_3d(Mpx(:,:,:,1),Mpy(:,:,:,1),Mpz(:,:,:,1),ng_f, &
                                      bpc(:,:,:,1),ng_c, &
                                      bp1(:,:,:,1),bp2(:,:,:,1),bp3(:,:,:,1),ng_e, &
                                      lo, hi)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine inverse_diag_lap

  subroutine inverse_diag_lap_2d(Minvx,Minvy,ng_f, &
                                 beta,ng_c, &
                                 beta_ed,ng_e,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_c,ng_e
    real(kind=dp_t), intent(inout) ::   Minvx(lo(1)-ng_f:,lo(2)-ng_f:)
    real(kind=dp_t), intent(inout) ::   Minvy(lo(1)-ng_f:,lo(2)-ng_f:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_c:,lo(2)-ng_c:)
    real(kind=dp_t), intent(in   ) :: beta_ed(lo(1)-ng_e:,lo(2)-ng_e:)
    
    ! local
    integer :: i,j
    real(kind=dp_t) :: b

    if (visc_type .eq. -1) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             Minvx(i,j) = 1.d0 / (beta(i,j)+beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))
          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             Minvy(i,j) = 1.d0 / (beta(i,j)+beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))
          end do
       end do

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2))

       Minvx = 1.d0 / (4.d0*b)
       Minvy = 1.d0 / (4.d0*b)

    else if (visc_type .eq. -2) then

       do j=lo(2),hi(2)
          do i=lo(1),hi(1)+1
             Minvx(i,j) = 1.d0 / (2.d0*beta(i,j)+2.d0*beta(i-1,j)+beta_ed(i,j+1)+beta_ed(i,j))
          end do
       end do

       do j=lo(2),hi(2)+1
          do i=lo(1),hi(1)
             Minvy(i,j) = 1.d0 / (2.d0*beta(i,j)+2.d0*beta(i,j-1)+beta_ed(i+1,j)+beta_ed(i,j))
          end do
       end do

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2))

       Minvx = 1.d0 / (6.d0*b)
       Minvy = 1.d0 / (6.d0*b)

    else if (abs(visc_type) .eq. 3) then

       call bl_error("inverse_diag_lap_2d does not support abs(visc_type) = 3")

    end if

  end subroutine inverse_diag_lap_2d

  subroutine inverse_diag_lap_3d(Minvx,Minvy,Minvz,ng_f, &
                                 beta,ng_c, &
                                 beta_xy,beta_xz,beta_yz,ng_e,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_f,ng_c,ng_e
    real(kind=dp_t), intent(inout) ::   Minvx(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    real(kind=dp_t), intent(inout) ::   Minvy(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    real(kind=dp_t), intent(inout) ::   Minvz(lo(1)-ng_f:,lo(2)-ng_f:,lo(3)-ng_f:)
    real(kind=dp_t), intent(in   ) ::    beta(lo(1)-ng_c:,lo(2)-ng_c:,lo(3)-ng_c:)
    real(kind=dp_t), intent(in   ) :: beta_xy(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_xz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)
    real(kind=dp_t), intent(in   ) :: beta_yz(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:)

    ! local
    integer :: i,j,k
    real(kind=dp_t) :: b

    if (visc_type .eq. -1) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                Minvx(i,j,k) = 1.d0 / &
                     ( beta(i,j,k)+beta(i-1,j,k) &
                      +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                      +beta_xz(i,j,k)+beta_xz(i,j,k+1) )
             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                Minvy(i,j,k) = 1.d0 / &
                     ( beta(i,j,k)+beta(i,j-1,k) &
                      +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                      +beta_yz(i,j,k)+beta_yz(i,j,k+1) )
             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                Minvz(i,j,k) = 1.d0 / &
                        ( beta(i,j,k)+beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )
                end do
             end do
          end do

    else if (visc_type .eq. 1) then

       b = beta(lo(1),lo(2),lo(3))

       Minvx = 1.d0 / (6.d0*b)
       Minvy = 1.d0 / (6.d0*b)
       Minvz = 1.d0 / (6.d0*b)

    else if (visc_type .eq. -2) then

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)+1
                Minvx(i,j,k) = 1.d0 / &
                     ( 2.d0*beta(i,j,k)+2.d0*beta(i-1,j,k) &
                      +beta_xy(i,j,k)+beta_xy(i,j+1,k) &
                      +beta_xz(i,j,k)+beta_xz(i,j,k+1) )
             end do
          end do
       end do

       do k=lo(3),hi(3)
          do j=lo(2),hi(2)+1
             do i=lo(1),hi(1)
                Minvy(i,j,k) = 1.d0 / &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j-1,k) &
                        +beta_xy(i,j,k)+beta_xy(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j,k+1) )
             end do
          end do
       end do

       do k=lo(3),hi(3)+1
          do j=lo(2),hi(2)
             do i=lo(1),hi(1)
                Minvz(i,j,k) = 1.d0 / &
                        ( 2.d0*beta(i,j,k)+2.d0*beta(i,j,k-1) &
                        +beta_xz(i,j,k)+beta_xz(i+1,j,k) &
                        +beta_yz(i,j,k)+beta_yz(i,j+1,k) )
                end do
             end do
          end do

    else if (visc_type .eq. 2) then

       b = beta(lo(1),lo(2),lo(3))

       Minvx = 1.d0 / (8.d0*b)
       Minvy = 1.d0 / (8.d0*b)
       Minvz = 1.d0 / (8.d0*b)

    else if (abs(visc_type) .eq. 3) then

       call bl_error("inverse_diag_lap_3d does not support abs(visc_type) = 3")

    end if

  end subroutine inverse_diag_lap_3d

end module inverse_diag_lap_module
