module compute_HSE_pres_module

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: n_cells, grav
  use multifab_physbc_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: compute_HSE_pres

contains

  ! initialize pressure using hydrostatic equilibrium
  ! assumes gravity is in the y-direction in 2D and 3D

  subroutine compute_HSE_pres(mla,sold,pold,dx,the_bc_tower)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: sold(:)
    type(multifab) , intent(inout) :: pold(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    real(kind=dp_t), allocatable :: rhoavg(:), p0(:)
    real(kind=dp_t), allocatable :: rhosum_proc(:), rhosum(:)
    integer :: ncell_wide, ncell_hi, ng_s, ng_p

    integer :: i,dm,n,nlevs,lo(mla%dim),hi(mla%dim)
    real(kind=dp_t), pointer :: sp(:,:,:,:)
    real(kind=dp_t), pointer :: pp(:,:,:,:)

    type(bl_prof_timer),save :: bpt
    
    call build(bpt,"compute_HSE_pres")

    nlevs = mla%nlevel
    dm = mla%dim

    ! assumes gravity is in y-direction in 2D and 3D
    ncell_hi = n_cells(2)

    allocate(     rhoavg(0:ncell_hi-1))
    allocate(         p0(0:ncell_hi-1))
    allocate(rhosum_proc(0:ncell_hi-1))
    allocate(rhosum     (0:ncell_hi-1))

    rhosum_proc = 0.d0

    if (dm .eq. 2) then
       ncell_wide = n_cells(1)
    else
       ncell_wide = n_cells(1)*n_cells(3)
    end if

    ng_s = sold(1)%ng
    ng_p = pold(1)%ng

    ! compute 1D array of average density as a function of height
    ! need to use parallel reduce tricks

    do n=1,nlevs
       do i=1,nfabs(sold(n))
          sp => dataptr(sold(n),i)
          lo = lwb(get_box(sold(n),i))
          hi = upb(get_box(sold(n),i))
          select case (dm)
          case (2)
             call sum_rho_2d(sp(:,:,1,1),rhosum_proc,lo,hi,ng_s)
          case (3)
             call sum_rho_3d(sp(:,:,:,1),rhosum_proc,lo,hi,ng_s)
          end select
       end do

       call parallel_reduce(rhosum(:),rhosum_proc(:),MPI_SUM)

       ! compute rhobar by normalizing rhosum
       do i=0,n_cells(2)-1
          rhoavg(i) = rhosum(i)/dble(ncell_wide)
       end do

    end do

    ! compute 1D array of pressure, assuming p=0 in first cell
    ! and integrating dp/dz = rhoavg*g
    p0(0) = 0.d0
    do i=1,ncell_hi-1
       p0(i) = p0(i-1) + 0.5d0*dx(1,2)*(rhoavg(i)+rhoavg(i-1))*grav(2)
    end do

    ! copy the 1D array into the full pressure multifab
    do n=1,nlevs
       do i=1,nfabs(pold(n))
          pp => dataptr(pold(n),i)
          lo = lwb(get_box(pold(n),i))
          hi = upb(get_box(pold(n),i))
          select case (dm)
          case (2)
             call copy_p0_2d(pp(:,:,1,1),p0,lo,hi,ng_p)
          case (3)
             call copy_p0_3d(pp(:,:,:,1),p0,lo,hi,ng_p)
          end select
       end do
    end do

    ! presure ghost cells
    do n=1,nlevs
       call multifab_fill_boundary(pold(n))
       call multifab_physbc(pold(n),1,pres_bc_comp,1,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    deallocate(rhoavg,p0)

    call destroy(bpt)

  contains
    
    subroutine sum_rho_2d(rho,rhosum,lo,hi,ng_s)

      integer         , intent(in   ) :: lo(:), hi(:), ng_s
      real (kind=dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:)
      real (kind=dp_t), intent(inout) :: rhosum(0:)
      
      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            rhosum(j) = rhosum(j) + rho(i,j)
         end do
      end do
      
    end subroutine sum_rho_2d
    
    subroutine sum_rho_3d(rho,rhosum,lo,hi,ng_s)

      integer         , intent(in   ) :: lo(:), hi(:), ng_s
      real (kind=dp_t), intent(in   ) :: rho(lo(1)-ng_s:,lo(2)-ng_s:,lo(3)-ng_s:)
      real (kind=dp_t), intent(inout) :: rhosum(0:)
      
      ! local
      integer :: i,j,k

      do j=lo(2),hi(2)
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            rhosum(j) = rhosum(j) + rho(i,j,k)
         end do
         end do
      end do
      
    end subroutine sum_rho_3d

    subroutine copy_p0_2d(pres,p0,lo,hi,ng_p)

      integer         , intent(in   ) :: lo(:), hi(:), ng_p
      real (kind=dp_t), intent(inout) :: pres(lo(1)-ng_p:,lo(2)-ng_p:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      
      ! local
      integer :: i,j

      do j=lo(2),hi(2)
         do i=lo(1),hi(1)
            pres(i,j) = p0(j)
         end do
      end do
      
    end subroutine copy_p0_2d

    subroutine copy_p0_3d(pres,p0,lo,hi,ng_p)

      integer         , intent(in   ) :: lo(:), hi(:), ng_p
      real (kind=dp_t), intent(inout) :: pres(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real (kind=dp_t), intent(in   ) :: p0(0:)
      
      ! local
      integer :: i,j,k

      do j=lo(2),hi(2)
         do k=lo(3),hi(3)
         do i=lo(1),hi(1)
            pres(i,j,k) = p0(j)
         end do
         end do
      end do
      
    end subroutine copy_p0_3d

  end subroutine compute_HSE_pres

end module compute_HSE_pres_module
