module mk_grav_force_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use zero_edgeval_module
  use probin_common_module, only: grav

  implicit none

  private

  public :: mk_grav_force

contains

  subroutine mk_grav_force(mla,m_force,s_fc_old,s_fc_new,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_force(:,:)
    type(multifab) , intent(in   ) :: s_fc_old(:,:)
    type(multifab) , intent(in   ) :: s_fc_new(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: i,n,ng_s,ng_u,dm,nlevs
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: fxp(:,:,:,:)
    real(kind=dp_t), pointer :: fyp(:,:,:,:)
    real(kind=dp_t), pointer :: fzp(:,:,:,:)
    real(kind=dp_t), pointer :: sox(:,:,:,:)
    real(kind=dp_t), pointer :: soy(:,:,:,:)
    real(kind=dp_t), pointer :: soz(:,:,:,:)
    real(kind=dp_t), pointer :: snx(:,:,:,:)
    real(kind=dp_t), pointer :: sny(:,:,:,:)
    real(kind=dp_t), pointer :: snz(:,:,:,:)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"mk_grav_force")

    nlevs = mla%nlevel
    dm = mla%dim

    ng_s  = s_fc_old(1,1)%ng
    ng_u  = m_force(1,1)%ng

    do n=1,nlevs
       do i=1, nfabs(m_force(n,1))
          fxp => dataptr(m_force(n,1), i)
          fyp => dataptr(m_force(n,2), i)
          sox => dataptr(s_fc_old(n,1), i)
          soy => dataptr(s_fc_old(n,2), i)
          snx => dataptr(s_fc_new(n,1), i)
          sny => dataptr(s_fc_new(n,2), i)
          lo = lwb(get_box(m_force(n,1), i))
          hi = upb(get_box(m_force(n,1), i))
          select case (dm)
          case (2)
             call mk_grav_force_2d(fxp(:,:,1,1), fyp(:,:,1,1), ng_u, &
                                   sox(:,:,1,:), soy(:,:,1,:), &
                                   snx(:,:,1,:), sny(:,:,1,:), ng_s, lo, hi)
          case (3)
             fzp => dataptr(m_force(n,3), i)
             soz => dataptr(s_fc_old(n,3), i)
             snz => dataptr(s_fc_new(n,3), i)
             call mk_grav_force_3d(fxp(:,:,:,1), fyp(:,:,:,1), fzp(:,:,:,1), ng_u, &
                                   sox(:,:,:,:), soy(:,:,:,:), soz(:,:,:,:), &
                                   snx(:,:,:,:), sny(:,:,:,:), snz(:,:,:,:), ng_s, lo, hi)
          end select
       end do

       ! zero wall boundary values
       call zero_edgeval_physical(m_force(n,:),1,1,the_bc_tower%bc_tower_array(n))

    enddo

    call destroy(bpt)

  end subroutine mk_grav_force

  subroutine mk_grav_force_2d(m_forcex,m_forcey,ng_u,rho_oldx,rho_oldy,rho_newx,rho_newy, &
                              ng_s,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rho_oldx(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_oldy(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_newx(lo(1)-ng_s:,lo(2)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_newy(lo(1)-ng_s:,lo(2)-ng_s:,:)

    ! local
    integer i,j

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
       m_forcex(i,j) = m_forcex(i,j) + &
            0.5d0*grav(1)*(rho_oldx(i,j,1)+rho_newx(i,j,1))
    end do
    end do

    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
       m_forcey(i,j) = m_forcey(i,j) + &
            0.5d0*grav(2)*(rho_oldy(i,j,1)+rho_newy(i,j,1))
    end do
    end do

  end subroutine mk_grav_force_2d

  subroutine mk_grav_force_3d(m_forcex,m_forcey,m_forcez,ng_u, &
                              rho_oldx,rho_oldy,rho_oldz, &
                              rho_newx,rho_newy,rho_newz,ng_s,lo,hi)

    integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_s
    real(kind=dp_t), intent(inout) :: m_forcex(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcey(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(inout) :: m_forcez(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
    real(kind=dp_t), intent(in   ) :: rho_oldx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_oldy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_oldz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_newx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_newy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)
    real(kind=dp_t), intent(in   ) :: rho_newz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:,:)

    ! local
    integer i,j,k

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)+1
       m_forcex(i,j,k) = m_forcex(i,j,k) + &
            0.5d0*grav(1)*(rho_oldx(i,j,k,1)+rho_newx(i,j,k,1))
    end do
    end do
    end do

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)+1
    do i=lo(1),hi(1)
       m_forcey(i,j,k) = m_forcey(i,j,k) + &
            0.5d0*grav(2)*(rho_oldy(i,j,k,1)+rho_newy(i,j,k,1))
    end do
    end do
    end do

    do k=lo(3),hi(3)+1
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       m_forcez(i,j,k) = m_forcez(i,j,k) + &
            0.5d0*grav(3)*(rho_oldz(i,j,k,1)+rho_newz(i,j,k,1))
    end do
    end do
    end do

  end subroutine mk_grav_force_3d

end module mk_grav_force_module
