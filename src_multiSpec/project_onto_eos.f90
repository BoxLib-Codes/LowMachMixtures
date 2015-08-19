module project_onto_eos_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use probin_common_module, only: rhobar
  use probin_gmres_module, only: mg_rel_tol
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: project_onto_eos

contains

  subroutine project_onto_eos(mla,rho)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)

    ! local
    integer i,n,nlevs,dm,ng_r
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: sp(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"project_onto_eos")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_r = rho(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          sp  => dataptr(rho(n), i)
          lo = lwb(get_box(rho(n), i))
          hi = upb(get_box(rho(n), i))
          select case (dm)
          case (2)
             call project_onto_eos_2d(sp(:,:,1,:), ng_r, lo, hi)
          case (3)
             call project_onto_eos_3d(sp(:,:,:,:), ng_r, lo, hi)
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine project_onto_eos

  subroutine project_onto_eos_2d(rho,ng_r,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,:)

    real(kind=dp_t) :: rho_tilde(lo(1):hi(1),lo(2):hi(2),nspecies)
    real(kind=dp_t) :: rhobar_sq,delta_eos,sum_spec(nspecies)
    real(kind=dp_t) :: rho_tmp,sum_change(nspecies),w(nspecies)
    integer i,j,l,ncell

    ! number of cells on the grid
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)

    sum_spec(:) = 0.d0
    sum_change(:) = 0.d0

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       ! compute mass fractions, w_i = rho_i/rho
       rho_tmp = 0.d0
       do l=1,nspecies
          rho_tmp = rho_tmp + rho(i,j,l)
       end do
       do l=1,nspecies
          w(l) = rho(i,j,l)/rho_tmp
       end do

       ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
       rhobar_sq = 0.d0
       do l=1,nspecies
          rhobar_sq = rhobar_sq + w(l)/rhobar(l)**2
       end do
       rhobar_sq = 1.d0/rhobar_sq
          
       ! delta_eos = sum_l (rho_l/rhobar_l) - 1
       delta_eos = -1.d0
       do l=1,nspecies
          delta_eos = delta_eos + rho(i,j,l)/rhobar(l)
       end do
       
       do l=1,nspecies
          ! rho_tilde_i = rho - w_i*(rhobar_sq/rhobar_i) * delta_eos
          rho_tilde(i,j,l) = rho(i,j,l) - w(l)*(rhobar_sq/rhobar(l))*delta_eos
          ! sum_spec_i = sum (rho_i - rho_tilde_i)
          sum_spec(l) = sum_spec(l) + rho(i,j,l) - rho_tilde(i,j,l)
       end do

    end do
    end do

    sum_spec(:) = sum_spec(:) / dble(ncell)

    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
    do l=1,nspecies
       rho_tmp = rho(i,j,l)
       rho(i,j,l) = rho_tilde(i,j,l) + sum_spec(l)
       sum_change(l) = sum_change(l) + (rho(i,j,l)-rho_tmp)**2
    end do
    end do
    end do

    ! redefine rhobar_sq = sum(rhobar_i^2)
    rhobar_sq = 0.d0
    do l=1,nspecies
       rhobar_sq = rhobar_sq + rhobar(l)**2
    end do
    sum_change(:) = sqrt(sum_change(:)/dble(ncell))
    sum_change(:) = sum_change(:) / sqrt(rhobar_sq)
    if(any( sum_change(1:nspecies) > 1000*mg_rel_tol)) then      
       call bl_warn('EOS adjustment exceeded Poisson solver tolerance')
       print*,sum_change(1:nspecies),mg_rel_tol
    end if

  end subroutine project_onto_eos_2d

  subroutine project_onto_eos_3d(rho,ng_r,lo,hi)

    ! L2 Projection onto EOS Constraint

    integer        , intent(in   ) :: lo(:),hi(:),ng_r
    real(kind=dp_t), intent(inout) :: rho(lo(1)-ng_r:,lo(2)-ng_r:,lo(3)-ng_r:,:)

    real(kind=dp_t) :: rho_tilde(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),nspecies)
    real(kind=dp_t) :: rhobar_sq,delta_eos,sum_spec(nspecies)
    real(kind=dp_t) :: rho_tmp,sum_change(nspecies),w(nspecies)
    integer i,j,k,l,ncell

    ! number of cells on the grid
    ncell = (hi(1)-lo(1)+1)*(hi(2)-lo(2)+1)*(hi(3)-lo(3)+1)

    sum_spec(:) = 0.d0
    sum_change(:) = 0.d0

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)

       ! compute mass fractions, w_i = rho_i/rho
       rho_tmp = 0.d0
       do l=1,nspecies
          rho_tmp = rho_tmp + rho(i,j,k,l)
       end do
       do l=1,nspecies
          w(l) = rho(i,j,k,l)/rho_tmp
       end do
          
       ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
       rhobar_sq = 0.d0
       do l=1,nspecies
          rhobar_sq = rhobar_sq + w(l)/rhobar(l)**2
       end do
       rhobar_sq = 1.d0/rhobar_sq

       ! delta_eos = sum_l (rho_l/rhobar_l) - 1
       delta_eos = -1.d0
       do l=1,nspecies
          delta_eos = delta_eos + rho(i,j,k,l)/rhobar(l)
       end do

       do l=1,nspecies
          ! rho_tilde_i = rho - w_i*(rhobar_sq/rhobar_i) * delta_eos
          rho_tilde(i,j,k,l) = rho(i,j,k,l) - w(l)*(rhobar_sq/rhobar(l))*delta_eos
          ! sum_spec_i = sum (rho_i - rho_tilde_i)
          sum_spec(l) = sum_spec(l) + rho(i,j,k,l) - rho_tilde(i,j,k,l)
       end do

    end do
    end do
    end do

    sum_spec(:) = sum_spec(:) / dble(ncell)

    do k=lo(3),hi(3)
    do j=lo(2),hi(2)
    do i=lo(1),hi(1)
       do l=1,nspecies
          rho_tmp = rho(i,j,k,l)
          rho(i,j,k,l) = rho_tilde(i,j,k,l) + sum_spec(l)
          sum_change(l) = sum_change(l) + (rho(i,j,k,l)-rho_tmp)**2
       end do
    end do
    end do
    end do

    ! redefine rhobar_sq = sum(rhobar_i^2)
    rhobar_sq = 0.d0
    do l=1,nspecies
       rhobar_sq = rhobar_sq + rhobar(l)**2
    end do
    sum_change(:) = sqrt(sum_change(:)/dble(ncell))
    sum_change(:) = sum_change(:) / sqrt(rhobar_sq)
    if(any( sum_change(1:nspecies) > 1000*mg_rel_tol)) then      
       call bl_warn('EOS adjustment exceeded Poisson solver tolerance')
       print*,sum_change(1:nspecies),mg_rel_tol
    end if


  end subroutine project_onto_eos_3d

end module project_onto_eos_module
