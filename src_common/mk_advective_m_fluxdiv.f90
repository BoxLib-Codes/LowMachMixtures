module mk_advective_m_fluxdiv_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_stag_module

  implicit none

  private

  public :: mk_advective_m_fluxdiv

contains

  subroutine mk_advective_m_fluxdiv(mla,umac,m,m_update,dx,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) ::     umac(:,:)
    type(multifab) , intent(in   ) ::        m(:,:)
    type(multifab) , intent(inout) :: m_update(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local
    integer :: i,n,nlevs,dm,ng_m,ng_u,ng_a
    integer :: lo(mla%dim),hi(mla%dim)

    real(kind=dp_t), pointer :: ump(:,:,:,:)
    real(kind=dp_t), pointer :: vmp(:,:,:,:)
    real(kind=dp_t), pointer :: wmp(:,:,:,:)
    real(kind=dp_t), pointer :: mxp(:,:,:,:)
    real(kind=dp_t), pointer :: myp(:,:,:,:)
    real(kind=dp_t), pointer :: mzp(:,:,:,:)
    real(kind=dp_t), pointer :: axp(:,:,:,:)
    real(kind=dp_t), pointer :: ayp(:,:,:,:)
    real(kind=dp_t), pointer :: azp(:,:,:,:)
    
    type(bl_prof_timer),save :: bpt

    call build(bpt,"mk_advective_m_fluxdiv")

    nlevs = mla%nlevel
    dm    = mla%dim

    ng_u = umac(1,1)%ng
    ng_m = m(1,1)%ng
    ng_a = m_update(1,1)%ng

    do n=1,nlevs
       do i=1,nfabs(m(n,1))
          ump => dataptr(umac(n,1), i)
          vmp => dataptr(umac(n,2), i)
          mxp => dataptr(m(n,1), i)
          myp => dataptr(m(n,2), i)
          axp => dataptr(m_update(n,1), i)
          ayp => dataptr(m_update(n,2), i)
          lo = lwb(get_box(m(n,1), i))
          hi = upb(get_box(m(n,1), i))
          select case (dm)
          case (2)
             call mk_advective_m_fluxdiv_2d(ump(:,:,1,1), vmp(:,:,1,1), ng_u, &
                                            mxp(:,:,1,1), myp(:,:,1,1), ng_m, &
                                            axp(:,:,1,1), ayp(:,:,1,1), ng_a, &
                                            lo, hi, dx(n,:))
          case (3)
             wmp => dataptr(umac(n,3), i)
             mzp => dataptr(m(n,3), i)
             azp => dataptr(m_update(n,3), i)
             call mk_advective_m_fluxdiv_3d(ump(:,:,:,1), vmp(:,:,:,1), wmp(:,:,:,1), ng_u, &
                                            mxp(:,:,:,1), myp(:,:,:,1), mzp(:,:,:,1), ng_m, &
                                            axp(:,:,:,1), ayp(:,:,:,1), azp(:,:,:,1), ng_a, &
                                            lo, hi, dx(n,:))
          end select
       end do

       do i=1,dm
          ! set advective update on physical domain boundaries to zero
          call multifab_physbc_domainvel(m_update(n,i),vel_bc_comp+i-1, &
                                         the_bc_level(n),dx(n,:))
       end do

    end do

    call destroy(bpt)

  contains

    subroutine mk_advective_m_fluxdiv_2d(umac,vmac,ng_u,mx,my,ng_m, &
                                         m_updatex,m_updatey,ng_a,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_m,ng_a
      real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_u:,lo(2)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        mx(lo(1)-ng_m:,lo(2)-ng_m:)
      real(kind=dp_t), intent(in   ) ::        my(lo(1)-ng_m:,lo(2)-ng_m:)
      real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_a:,lo(2)-ng_a:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j

      real(kind=dp_t) :: mx_fluxx  (lo(1):hi(1)+2,lo(2):hi(2)  )
      real(kind=dp_t) :: mx_fluxy  (lo(1):hi(1)+1,lo(2):hi(2)+1)
      real(kind=dp_t) :: mx_fluxdiv(lo(1):hi(1)+1,lo(2):hi(2)  )

      real(kind=dp_t) ::   my_fluxx(lo(1):hi(1)+1,lo(2):hi(2)+1)
      real(kind=dp_t) ::   my_fluxy(lo(1):hi(1)  ,lo(2):hi(2)+2)
      real(kind=dp_t) :: my_fluxdiv(lo(1):hi(1)  ,lo(2):hi(2)+1)

      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      !=============================
      ! mx fluxes and divergence
      !=============================
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+2
            mx_fluxx(i,j) = 0.25d0*(mx(i-1,j)+mx(i,j))*(umac(i-1,j)+umac(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            mx_fluxy(i,j) = 0.25d0*(mx(i,j-1)+mx(i,j))*(vmac(i-1,j)+vmac(i,j))
         end do
      end do

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            mx_fluxdiv(i,j) = -( (mx_fluxx(i+1,j)-mx_fluxx(i,j)) * dxinv + &
                                 (mx_fluxy(i,j+1)-mx_fluxy(i,j)) * dxinv )
         end do
      end do

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
            m_updatex(i,j) = m_updatex(i,j) + mx_fluxdiv(i,j)
         end do
      end do

      !=============================
      ! my fluxes and divergence
      !=============================
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)+1
            my_fluxx(i,j) = 0.25d0*(my(i-1,j)+my(i,j))*(umac(i,j-1)+umac(i,j))
         end do
      end do

      do j=lo(2),hi(2)+2
         do i=lo(1),hi(1)
            my_fluxy(i,j) = 0.25d0*(my(i,j-1)+my(i,j))*(vmac(i,j-1)+vmac(i,j))
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            my_fluxdiv(i,j) = -( (my_fluxx(i+1,j)-my_fluxx(i,j)) * dxinv + &
                                 (my_fluxy(i,j+1)-my_fluxy(i,j)) * dxinv )
         end do
      end do

      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)
            m_updatey(i,j) = m_updatey(i,j) + my_fluxdiv(i,j)
         end do
      end do

    end subroutine mk_advective_m_fluxdiv_2d

    subroutine mk_advective_m_fluxdiv_3d(umac,vmac,wmac,ng_u,mx,my,mz,ng_m, &
                                         m_updatex,m_updatey,m_updatez,ng_a,lo,hi,dx)

      integer        , intent(in   ) :: lo(:),hi(:),ng_u,ng_m,ng_a
      real(kind=dp_t), intent(in   ) ::      umac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::      vmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::      wmac(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:)
      real(kind=dp_t), intent(in   ) ::        mx(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(in   ) ::        my(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(in   ) ::        mz(lo(1)-ng_m:,lo(2)-ng_m:,lo(3)-ng_m:)
      real(kind=dp_t), intent(inout) :: m_updatex(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout) :: m_updatey(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(inout) :: m_updatez(lo(1)-ng_a:,lo(2)-ng_a:,lo(3)-ng_a:)
      real(kind=dp_t), intent(in   ) :: dx(:)

      ! local
      integer :: i,j,k

      real(kind=dp_t) ::   mx_fluxx(lo(1):hi(1)+2,lo(2):hi(2)  ,lo(3):hi(3)  )
      real(kind=dp_t) ::   mx_fluxy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) ::   mx_fluxz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) :: mx_fluxdiv(lo(1):hi(1)+1,lo(2):hi(2)  ,lo(3):hi(3)  )

      real(kind=dp_t) ::   my_fluxx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) ::   my_fluxy(lo(1):hi(1)  ,lo(2):hi(2)+2,lo(3):hi(3)  )
      real(kind=dp_t) ::   my_fluxz(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) :: my_fluxdiv(lo(1):hi(1)  ,lo(2):hi(2)+1,lo(3):hi(3)  )

      real(kind=dp_t) ::   mz_fluxx(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) ::   mz_fluxy(lo(1):hi(1)+1,lo(2):hi(2)+1,lo(3):hi(3)+1)
      real(kind=dp_t) ::   mz_fluxz(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+2)
      real(kind=dp_t) :: mz_fluxdiv(lo(1):hi(1)  ,lo(2):hi(2)  ,lo(3):hi(3)+1)

      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      !=============================
      ! mx fluxes and divergence
      !=============================
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+2
               mx_fluxx(i,j,k) = 0.25d0*(mx(i-1,j,k)+mx(i,j,k))*(umac(i-1,j,k)+umac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               mx_fluxy(i,j,k) = 0.25d0*(mx(i,j-1,k)+mx(i,j,k))*(vmac(i-1,j,k)+vmac(i,j,k))
               mx_fluxz(i,j,k) = 0.25d0*(mx(i,j,k-1)+mx(i,j,k))*(wmac(i-1,j,k)+wmac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               mx_fluxdiv(i,j,k) = -( (mx_fluxx(i+1,j,k)-mx_fluxx(i,j,k)) * dxinv + &
                                      (mx_fluxy(i,j+1,k)-mx_fluxy(i,j,k)) * dxinv + &
                                      (mx_fluxz(i,j,k+1)-mx_fluxz(i,j,k)) * dxinv )
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
               m_updatex(i,j,k) = m_updatex(i,j,k) + mx_fluxdiv(i,j,k)
            end do
         end do
      end do

      !=============================
      ! my fluxes and divergence
      !=============================
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+2
            do i=lo(1),hi(1)
               my_fluxy(i,j,k) = 0.25d0*(my(i,j-1,k)+my(i,j,k))*(vmac(i,j-1,k)+vmac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               my_fluxx(i,j,k) = 0.25d0*(my(i-1,j,k)+my(i,j,k))*(umac(i,j-1,k)+umac(i,j,k))
               my_fluxz(i,j,k) = 0.25d0*(my(i,j,k-1)+my(i,j,k))*(wmac(i,j-1,k)+wmac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               my_fluxdiv(i,j,k) = -( (my_fluxx(i+1,j,k)-my_fluxx(i,j,k)) * dxinv + &
                                      (my_fluxy(i,j+1,k)-my_fluxy(i,j,k)) * dxinv + &
                                      (my_fluxz(i,j,k+1)-my_fluxz(i,j,k)) * dxinv )
            end do
         end do
      end do

      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
               m_updatey(i,j,k) = m_updatey(i,j,k) + my_fluxdiv(i,j,k)
            end do
         end do
      end do

      !=============================
      ! mz fluxes and divergence
      !=============================
      do k=lo(3),hi(3)+2
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               mz_fluxz(i,j,k) = 0.25d0*(mz(i,j,k-1)+mz(i,j,k))*(wmac(i,j,k-1)+wmac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)+1
               mz_fluxx(i,j,k) = 0.25d0*(mz(i-1,j,k)+mz(i,j,k))*(umac(i,j,k-1)+umac(i,j,k))
               mz_fluxy(i,j,k) = 0.25d0*(mz(i,j-1,k)+mz(i,j,k))*(vmac(i,j,k-1)+vmac(i,j,k))
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               mz_fluxdiv(i,j,k) = -( (mz_fluxx(i+1,j,k)-mz_fluxx(i,j,k)) * dxinv + &
                                      (mz_fluxy(i,j+1,k)-mz_fluxy(i,j,k)) * dxinv + &
                                      (mz_fluxz(i,j,k+1)-mz_fluxz(i,j,k)) * dxinv )
            end do
         end do
      end do

      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
               m_updatez(i,j,k) = m_updatez(i,j,k) + mz_fluxdiv(i,j,k)
            end do
         end do
      end do

    end subroutine mk_advective_m_fluxdiv_3d
    
  end subroutine mk_advective_m_fluxdiv

end module mk_advective_m_fluxdiv_module
