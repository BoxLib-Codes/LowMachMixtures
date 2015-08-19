module div_and_grad_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module

  implicit none

  private

  public :: compute_grad, compute_div

contains

  subroutine compute_grad(mla,phi,gradp,dx, &
                          start_incomp,start_bccomp,start_outcomp,num_comp, &
                          the_bc_level,increment_bccomp_in)

    ! compute the face-centered gradient of a cell-centered field

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi(:)
    type(multifab) , intent(inout) :: gradp(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: start_incomp,start_bccomp,start_outcomp,num_comp
    type(bc_level) , intent(in   ) :: the_bc_level(:)
    logical,  intent(in), optional :: increment_bccomp_in

    ! local
    integer :: n,i,dm,nlevs,ng_p,ng_g,comp,bccomp,outcomp
    integer :: lo(mla%dim),hi(mla%dim)
    logical :: increment_bccomp

    real(kind=dp_t), pointer :: pp(:,:,:,:)
    real(kind=dp_t), pointer :: gpx(:,:,:,:)
    real(kind=dp_t), pointer :: gpy(:,:,:,:)
    real(kind=dp_t), pointer :: gpz(:,:,:,:)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"compute_grad")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_p = phi(1)%ng
    ng_g = gradp(1,1)%ng

    increment_bccomp = .true.
    if (present(increment_bccomp_in)) then
       increment_bccomp = increment_bccomp_in
    end if
    
    !$omp parallel private(mfi,n,i,xnodalbox,ynodalbox,znodalbox,xlo,ylo,zlo) &
    !$omp private(xhi,yhi,zhi,pp,gpx,gpy,gpz,lo,hi,comp,bccomp,outcomp)
    do n=1,nlevs
       call mfiter_build(mfi, phi(n), tiling=.true.)

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

 !      do i=1,nfabs(phi(n))
          pp  => dataptr(phi(n), i)
          gpx => dataptr(gradp(n,1), i)
          gpy => dataptr(gradp(n,2), i)
          lo = lwb(get_box(phi(n), i))
          hi = upb(get_box(phi(n), i))
          
          do comp=start_incomp,start_incomp+num_comp-1
             if (increment_bccomp) then
                bccomp = start_bccomp + (comp-start_incomp)
             else
                bccomp = start_bccomp
             end if
             outcomp = start_outcomp + (comp-start_incomp)
             select case (dm)
             case (2)
                call compute_grad_2d(pp(:,:,1,comp), ng_p, &
                                     gpx(:,:,1,outcomp), gpy(:,:,1,outcomp), ng_g, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp),xlo,xhi,ylo,yhi)
             case (3)
                gpz => dataptr(gradp(n,3), i)
                call compute_grad_3d(pp(:,:,:,comp), ng_p, &
                                     gpx(:,:,:,outcomp), gpy(:,:,:,outcomp), gpz(:,:,:,outcomp), ng_g, &
                                     lo, hi, dx(n,:), &
                                     the_bc_level(n)%adv_bc_level_array(i,:,:,bccomp),xlo,xhi,ylo,yhi,zlo,zhi)
             end select
          end do
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  contains
    
    subroutine compute_grad_2d(phi,ng_p,gpx,gpy,ng_g,lo,hi,dx,bc,xlo,xhi,ylo,yhi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
      real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) :: gpx(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpy(lo(1)-ng_g:,lo(2)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      ! local
      integer :: i,j
      real(kind=dp_t) :: dxinv,twodxinv

      dxinv = 1.d0/dx(1)
      twodxinv = 2.d0*dxinv

      ! x-faces
      do j=xlo(2),xhi(2)
         do i=xlo(1),xhi(1)
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) * dxinv
         end do
      end do
   
      ! alter stencil at boundary since ghost value represents value at boundary
      if (xlo(1) .eq. lo(1)) then
      if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
         i=lo(1)
         do j=lo(2),hi(2)
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) * twodxinv
         end do
      end if
      end if
      if (xhi(1) .eq. hi(1)+1) then
      if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
         i=hi(1)+1
         do j=lo(2),hi(2)
            gpx(i,j) = ( phi(i,j)-phi(i-1,j) ) * twodxinv
         end do
      end if
      end if

      ! y-faces
      do j=ylo(2),yhi(2)
         do i=ylo(1),yhi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) * dxinv
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (ylo(2) .eq. lo(2)) then
      if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
         j=lo(2)
         do i=lo(1),hi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) * twodxinv
         end do
      end if
      end if

      ! alter stencil at boundary since ghost value represents value at boundary
      if (yhi(2) .eq. hi(2)+1) then
      if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
         j=hi(2)+1
         do i=lo(1),hi(1)
            gpy(i,j) = ( phi(i,j)-phi(i,j-1) ) * twodxinv
         end do
      end if
      end if

    end subroutine compute_grad_2d

    subroutine compute_grad_3d(phi,ng_p,gpx,gpy,gpz,ng_g,lo,hi,dx,bc,xlo,xhi,ylo,yhi,zlo,zhi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      integer        , intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
      real(kind=dp_t), intent(in   ) :: phi(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) :: gpx(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpy(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(inout) :: gpz(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      integer        , intent(in   ) :: bc(:,:)

      ! local
      integer :: i,j,k
      real(kind=dp_t) :: dxinv,twodxinv
      
      dxinv = 1.d0/dx(1)
      twodxinv = 2.d0*dxinv


      ! x-faces
      do k=xlo(3),xhi(3)
         do j=xlo(2),xhi(2)
            do i=xlo(1),xhi(1)
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) * dxinv
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (xlo(1) .eq. lo(1)) then 
      if (bc(1,1) .eq. FOEXTRAP .or. bc(1,1) .eq. HOEXTRAP .or. bc(1,1) .eq. EXT_DIR) then
         i=xlo(1)
         do k=xlo(3),xhi(3)
            do j=xlo(2),xhi(2)
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) * twodxinv
            end do
         end do
      end if
      end if
      if (xhi(1) .eq. hi(1)+1) then
      if (bc(1,2) .eq. FOEXTRAP .or. bc(1,2) .eq. HOEXTRAP .or. bc(1,2) .eq. EXT_DIR) then
         i=xhi(1)
         do k=xlo(3),xhi(3)
            do j=xlo(2),xhi(2)
               gpx(i,j,k) = ( phi(i,j,k)-phi(i-1,j,k) ) * twodxinv
            end do
         end do
      end if
      end if

      ! y-faces
      do k=ylo(3),yhi(3)
         do j=ylo(2),yhi(2)
            do i=ylo(1),yhi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) * dxinv
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (ylo(2) .eq. lo(2)) then
      if (bc(2,1) .eq. FOEXTRAP .or. bc(2,1) .eq. HOEXTRAP .or. bc(2,1) .eq. EXT_DIR) then
         j=ylo(2)
         do k=ylo(3),yhi(3)
            do i=ylo(1),yhi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) * twodxinv
            end do
         end do
      end if
      end if
      if (yhi(2) .eq. hi(2)+1) then
      if (bc(2,2) .eq. FOEXTRAP .or. bc(2,2) .eq. HOEXTRAP .or. bc(2,2) .eq. EXT_DIR) then
         j=yhi(2)
         do k=ylo(3),yhi(3)
            do i=ylo(1),yhi(1)
               gpy(i,j,k) = ( phi(i,j,k)-phi(i,j-1,k) ) * twodxinv
            end do
         end do
      end if
      end if

      ! z-faces
      do k=zlo(3),zhi(3)
         do j=zlo(2),zhi(2)
            do i=zlo(1),zhi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) * dxinv
            end do
         end do
      end do

      ! alter stencil at boundary since ghost value represents value at boundary
      if (zlo(3) .eq. lo(3)) then
      if (bc(3,1) .eq. FOEXTRAP .or. bc(3,1) .eq. HOEXTRAP .or. bc(3,1) .eq. EXT_DIR) then
         k=zlo(3)
         do j=zlo(2),zhi(2)
            do i=zlo(1),zhi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) * twodxinv
            end do
         end do
      end if
      end if

      ! alter stencil at boundary since ghost value represents value at boundary
      if (zhi(3) .eq. hi(3)+1) then
      if (bc(3,2) .eq. FOEXTRAP .or. bc(3,2) .eq. HOEXTRAP .or. bc(3,2) .eq. EXT_DIR) then
         k=zhi(3)
         do j=zlo(2),zhi(2)
            do i=zlo(1),zhi(1)
               gpz(i,j,k) = ( phi(i,j,k)-phi(i,j,k-1) ) * twodxinv
            end do
         end do
      end if
      end if

    end subroutine compute_grad_3d

  end subroutine compute_grad

  subroutine compute_div(mla,phi_fc,div,dx,start_incomp,start_outcomp,num_comp,increment_in)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: phi_fc(:,:)
    type(multifab) , intent(inout) :: div(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    integer        , intent(in   ) :: start_incomp, start_outcomp, num_comp
    logical,  intent(in), optional :: increment_in

    real(kind=dp_t), pointer :: pxp(:,:,:,:) 
    real(kind=dp_t), pointer :: pyp(:,:,:,:) 
    real(kind=dp_t), pointer :: pzp(:,:,:,:) 
    real(kind=dp_t), pointer :: dp(:,:,:,:) 
    integer :: i,n,nlevs,dm,ng_p,ng_d,comp,outcomp
    integer :: lo(mla%dim),hi(mla%dim)
    logical :: increment

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"compute_div")

    ! do we increment or overwrite div?
    increment = .false.
    if (present(increment_in)) increment = increment_in

    dm = mla%dim
    nlevs = mla%nlevel

    ng_p = phi_fc(1,1)%ng
    ng_d = div(1)%ng

    !$omp parallel private(mfi,n,i,tilebox,tlo,thi,pxp,pyp,pzp,dp,lo,hi,comp,outcomp)
    do n = 1,nlevs
       call mfiter_build(mfi, div(n), tiling=.true.)

     do while (more_tile(mfi))
       i = get_fab_index(mfi)

       tilebox = get_tilebox(mfi)
       tlo = lwb(tilebox)
       thi = upb(tilebox)

!       do i = 1, nfabs(div(n))
          pxp => dataptr(phi_fc(n,1), i)
          pyp => dataptr(phi_fc(n,2), i)
          dp => dataptr(div(n), i)
          lo =  lwb(get_box(div(n), i))
          hi =  upb(get_box(div(n), i))
          do comp=start_incomp,start_incomp+num_comp-1
             outcomp = start_outcomp + (comp-start_incomp)
             select case (dm)
             case (2)
                call compute_div_2d(pxp(:,:,1,comp), pyp(:,:,1,comp), ng_p, &
                                    dp(:,:,1,outcomp), ng_d, dx(n,:),lo, hi, increment,tlo,thi)
             case (3) 
                pzp => dataptr(phi_fc(n,3), i)
                call compute_div_3d(pxp(:,:,:,comp), pyp(:,:,:,comp), pzp(:,:,:,comp), ng_p, &
                                    dp(:,:,:,outcomp), ng_d, dx(n,:),lo, hi, increment,tlo,thi)
             end select
          end do
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  contains

    subroutine compute_div_2d(phix,phiy,ng_p,div,ng_d,dx,lo,hi,increment,tlo,thi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_d,tlo(:),thi(:)
      real(kind=dp_t), intent(in   ) :: phix(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiy(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) ::  div(lo(1)-ng_d:,lo(2)-ng_d:)
      real(kind=dp_t), intent(in   ) ::   dx(:)
      logical        , intent(in   ) :: increment

      integer :: i,j
      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      if (increment) then

         do j = tlo(2),thi(2)
         do i = tlo(1),thi(1)
            div(i,j) = div(i,j) + &
                 (phix(i+1,j) - phix(i,j)) * dxinv  + &
                 (phiy(i,j+1) - phiy(i,j)) * dxinv
         end do
         end do

      else

         do j = tlo(2),thi(2)
         do i = tlo(1),thi(1)
            div(i,j) = &
                 (phix(i+1,j) - phix(i,j)) * dxinv + &
                 (phiy(i,j+1) - phiy(i,j)) * dxinv
         end do
         end do
         
      end if

    end subroutine compute_div_2d

    subroutine compute_div_3d(phix,phiy,phiz,ng_p,div,ng_d,dx,lo,hi,increment,tlo,thi)

      integer        , intent(in   ) :: lo(:),hi(:),ng_p,ng_d,tlo(:),thi(:)
      real(kind=dp_t), intent(in   ) :: phix(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiy(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(in   ) :: phiz(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) ::  div(lo(1)-ng_d:,lo(2)-ng_d:,lo(3)-ng_d:)
      real(kind=dp_t), intent(in   ) :: dx(:)
      logical        , intent(in   ) :: increment

      integer :: i,j,k
      real(kind=dp_t) :: dxinv

      dxinv = 1.d0/dx(1)

      if (increment) then
         do k = tlo(3),thi(3)
         do j = tlo(2),thi(2)
         do i = tlo(1),thi(1)
            div(i,j,k) = div(i,j,k) + &
                 (phix(i+1,j,k) - phix(i,j,k)) * dxinv + &
                 (phiy(i,j+1,k) - phiy(i,j,k)) * dxinv + &
                 (phiz(i,j,k+1) - phiz(i,j,k)) * dxinv
         end do
         end do
         end do
      else
         do k = tlo(3),thi(3)
         do j = tlo(2),thi(2)
         do i = tlo(1),thi(1)
            div(i,j,k) = &
                 (phix(i+1,j,k) - phix(i,j,k)) * dxinv + &
                 (phiy(i,j+1,k) - phiy(i,j,k)) * dxinv + &
                 (phiz(i,j,k+1) - phiz(i,j,k)) * dxinv
         end do
         end do
         end do
      end if

    end subroutine compute_div_3d

  end subroutine compute_div

end module div_and_grad_module
