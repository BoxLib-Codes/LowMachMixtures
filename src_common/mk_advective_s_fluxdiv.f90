module mk_advective_s_fluxdiv_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_stag_module

  implicit none

  private

  public :: mk_advective_s_fluxdiv

contains

  subroutine mk_advective_s_fluxdiv(mla,umac,s_fc,s_update,dx,start_comp,num_comp)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: s_fc(:,:)
    type(multifab) , intent(inout) :: s_update(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: start_comp,num_comp

    ! local
    integer :: i,n,nlevs,dm,ng_s,ng_u,ng_a
    integer :: lo(mla%dim),hi(mla%dim),comp

    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: ap(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    
    type(bl_prof_timer),save :: bpt

    call build(bpt,"mk_advective_s_fluxdiv")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_s = s_fc(1,1)%ng
    ng_u = umac(1,1)%ng
    ng_a = s_update(1)%ng

    do n=1,nlevs
       do i=1,nfabs(s_update(n))
          spx => dataptr(s_fc(n,1), i)
          spy => dataptr(s_fc(n,2), i)
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          ap  => dataptr(s_update(n), i)
          lo = lwb(get_box(s_update(n), i))
          hi = upb(get_box(s_update(n), i))
          do comp=start_comp,start_comp+num_comp-1
             select case (dm)
             case (2)
                call mk_advective_s_fluxdiv_2d(spx(:,:,1,comp), spy(:,:,1,comp), &
                                               ump(:,:,1,1), vmp(:,:,1,1), &
                                               ap(:,:,1,comp), ng_s, ng_u, ng_a, lo, hi, dx(n,:))
             case (3)
                wmp => dataptr(umac(n,3), i)
                spz => dataptr(s_fc(n,3), i)
                call mk_advective_s_fluxdiv_3d(spx(:,:,:,comp), spy(:,:,:,comp), spz(:,:,:,comp), &
                                               ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), &
                                               ap(:,:,:,comp), ng_s, ng_u, ng_a, lo, hi, dx(n,:))
             end select
          end do
       end do
    end do

    call destroy(bpt)

  contains

    subroutine mk_advective_s_fluxdiv_2d(sx,sy,umac,vmac,s_update,ng_s,ng_u,ng_a,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u,ng_a
      real(kind=dp_t), intent(in   ) ::       sx(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::       sy(lo(1)-ng_s:,lo(2)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     umac(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::     vmac(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1)

      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            fluxx(i,j) = umac(i,j)*sx(i,j)
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            fluxy(i,j) = vmac(i,j)*sy(i,j)
         end do
      end do

      !=============================
      ! Calculate the divergence of the advective flux:
      !=============================
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            s_update(i,j) = s_update(i,j) - ( &
                 (fluxx(i+1,j)-fluxx(i,j)) * dxinv &
                 + (fluxy(i,j+1)-fluxy(i,j)) * dxinv )
         end do
      end do

    end subroutine mk_advective_s_fluxdiv_2d

    subroutine mk_advective_s_fluxdiv_3d(sx,sy,sz,umac,vmac,wmac,s_update,ng_s,ng_u,ng_a,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_s,ng_u,ng_a
      real(kind=dp_t), intent(in   ) ::       sx(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::       sy(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::       sz(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real(kind=dp_t), intent(in   ) ::     umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::     vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::     wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(inout) :: s_update(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k

      real(kind=dp_t) :: fluxx(lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
      real(kind=dp_t) :: fluxy(lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
      real(kind=dp_t) :: fluxz(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)

      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)
      
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               fluxx(i,j,k) = umac(i,j,k)*sx(i,j,k)
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               fluxy(i,j,k) = vmac(i,j,k)*sy(i,j,k)
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               fluxz(i,j,k) = wmac(i,j,k)*sz(i,j,k)
            end do
         end do
      end do

      !=============================
      ! Calculate the divergence of the advective flux:
      !=============================
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               s_update(i,j,k) = s_update(i,j,k) - ( &
                    (fluxx(i+1,j,k)-fluxx(i,j,k)) * dxinv &
                    + (fluxy(i,j+1,k)-fluxy(i,j,k)) * dxinv &
                    + (fluxz(i,j,k+1)-fluxz(i,j,k)) * dxinv )
            end do
         end do
      end do

    end subroutine mk_advective_s_fluxdiv_3d

  end subroutine mk_advective_s_fluxdiv

end module mk_advective_s_fluxdiv_module
