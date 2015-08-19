module correction_flux_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use probin_multispecies_module, only: fraction_tolerance, nspecies

  implicit none

  private

  public :: correction_flux

contains

  subroutine correction_flux(mla,rho,rhotot,flux,the_bc_level)

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: flux(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: n,i,dm,nlevs,ng_p,ng_g
    integer :: lo(mla%dim),hi(mla%dim)

    ! pointer for rho, rhotot, flux_x, flux_y and flux_z 
    real(kind=dp_t), pointer :: dp(:,:,:,:)      ! for rho
    real(kind=dp_t), pointer :: dp1(:,:,:,:)     ! for rhotot
    real(kind=dp_t), pointer :: flux_x(:,:,:,:)
    real(kind=dp_t), pointer :: flux_y(:,:,:,:)
    real(kind=dp_t), pointer :: flux_z(:,:,:,:)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "correction_flux")

    dm    = mla%dim       ! dimensionality
    nlevs = mla%nlevel    ! number of levels 
    ng_p  = rho(1)%ng     ! number of ghost cells for rho
    ng_g  = flux(1,1)%ng  ! number of ghost cells for flux

    ! loop over all boxes   
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp  => dataptr(rho(n), i)
          dp1 => dataptr(rhotot(n), i)
          flux_x => dataptr(flux(n,1), i)
          flux_y => dataptr(flux(n,2), i)
          lo = lwb(get_box(rho(n), i))
          hi = upb(get_box(rho(n), i))
          
          select case (dm)
          case (2)
          call correction_flux_2d(dp(:,:,1,:),dp1(:,:,1,1),ng_p,flux_x(:,:,1,:),flux_y(:,:,1,:),&
                                  ng_g,lo,hi)
          case (3)
          flux_z => dataptr(flux(n,3), i)
          call correction_flux_3d(dp(:,:,:,:),dp1(:,:,:,1),ng_p,flux_x(:,:,:,:),flux_y(:,:,:,:),&
                                  flux_z(:,:,:,:),ng_g,lo,hi)
          end select
       end do
    end do

    call destroy(bpt)

  contains
    
    subroutine correction_flux_2d(rho,rhotot,ng_p,flux_x,flux_y,ng_g,lo,hi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_p:,lo(2)-ng_p:,:)
      real(kind=dp_t), intent(in   ) :: rhotot(lo(1)-ng_p:,lo(2)-ng_p:)
      real(kind=dp_t), intent(inout) :: flux_x(lo(1)-ng_g:,lo(2)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_y(lo(1)-ng_g:,lo(2)-ng_g:,:)

      ! local
      integer         :: i,j,n
      real(kind=dp_t) :: sumx,sumy,corr,total_corr

      ! x-faces
      total_corr=0.0d0
      do j=lo(2),hi(2)
         do i=lo(1),hi(1)+1
               
            ! free the data
            sumx = 0.d0
            corr = 0.d0
              
            ! sum the x-fluxes upto nspecies-1 
            do n=1, nspecies-1
               sumx = sumx + flux_x(i,j,n)
            end do
              
            ! caculate corr and print error if not zero 
            corr = flux_x(i,j,nspecies) + sumx
            
            if(corr .gt. rhotot(i,j)*fraction_tolerance) then
               !write(*,*) "Error: sum of x-fluxes = ", corr, " fluxes=", flux_x(i,j,:)
            end if
              
            ! correct x-flux for last species  
            flux_x(i,j,nspecies) = -sumx     
            total_corr = total_corr + abs(corr)
            !if(i .eq. 32 .and. j.gt.13 .and. j.lt.19) print*, 'x-flux=', flux_x(i,j,:)

         end do
      end do
      !write(*,*) "x flux correction = ", total_corr
   
      ! y-faces
      total_corr=0.0d0
      do j=lo(2),hi(2)+1
         do i=lo(1),hi(1)

            ! free the data
            sumy  = 0.d0
            corr = 0.d0
              
            ! sum the y-fluxes upto nspecies-1 
            do n=1, nspecies-1
               sumy = sumy + flux_y(i,j,n)
            end do
              
            ! caculate corr and print error if not zero 
            corr = flux_y(i,j,nspecies) + sumy
            
            if(corr .gt. rhotot(i,j)*fraction_tolerance) then
               !write(*,*) "Error: sum of y-fluxes = ", corr, " fluxes=", flux_y(i,j,:)
            end if
              
            ! correct y-flux for last species  
            flux_y(i,j,nspecies) = -sumy             
            total_corr = total_corr + abs(corr)
            !if(i .eq. 32 .and. j.gt.14 .and. j.lt.18) print*, 'y-flux=', flux_y(i,j,:)

         end do
      end do
     !write(*,*) "y flux correction = ", total_corr

    end subroutine correction_flux_2d

    subroutine correction_flux_3d(rho,rhotot,ng_p,flux_x,flux_y,flux_z,ng_g,lo,hi)

      integer        , intent(in   ) :: ng_p,ng_g,lo(:),hi(:)
      real(kind=dp_t), intent(in   ) :: rho(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:,:)
      real(kind=dp_t), intent(in   ) :: rhotot(lo(1)-ng_p:,lo(2)-ng_p:,lo(3)-ng_p:)
      real(kind=dp_t), intent(inout) :: flux_x(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_y(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)
      real(kind=dp_t), intent(inout) :: flux_z(lo(1)-ng_g:,lo(2)-ng_g:,lo(3)-ng_g:,:)

      ! local
      integer         :: i,j,k,n
      real(kind=dp_t) :: sumx,sumy,sumz,corr,total_corr
      
      ! x-faces
      total_corr=0.0d0
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)+1
   
               ! free the data
               sumx = 0.d0
               corr = 0.d0
             
               ! sum the x-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumx = sumx + flux_x(i,j,k,n)
               end do
              
               ! caculate corr and print error if not zero 
               corr = flux_x(i,j,k,nspecies) + sumx 
              
               if(corr .gt. rhotot(i,j,k)*fraction_tolerance) then
                  !write(*,*) "Error: sum of x-fluxes = ", corr, " fluxes=", flux_x(i,j,k,:)
               end if
              
               ! correct x-flux for last species  
               flux_x(i,j,k,nspecies) = -sumx             
               total_corr = total_corr + abs(corr) 

            end do
         end do
      end do
      !write(*,*) "x flux correction = ", total_corr
      
      ! y-faces
      total_corr=0.0d0
      do k=lo(3),hi(3)
         do j=lo(2),hi(2)+1
            do i=lo(1),hi(1)
              
               ! free the data
               sumy  = 0.d0
               corr = 0.d0
              
               ! sum the y-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumy = sumy + flux_y(i,j,k,n)
               end do
              
               ! caculate corr and print error if not zero 
               corr = flux_y(i,j,k,nspecies) + sumy 
               
               if(corr .gt. rhotot(i,j,k)*fraction_tolerance) then
                  !write(*,*) "Error: sum of y-fluxes = ", corr, " fluxes=", flux_y(i,j,k,:)
               end if
              
               ! correct y-flux for last species  
               flux_y(i,j,k,nspecies) = -sumy             
               total_corr = total_corr + abs(corr) 
 
            end do
         end do
      end do
      !write(*,*) "y flux correction = ", total_corr

      ! z-faces
      total_corr=0.0d0
      do k=lo(3),hi(3)+1
         do j=lo(2),hi(2)
            do i=lo(1),hi(1)
            
               ! free the data
               sumz  = 0.d0
               corr = 0.d0
              
               ! sum the z-fluxes upto nspecies-1 
               do n=1, nspecies-1
                  sumz = sumz + flux_z(i,j,k,n)
               end do
              
               ! caculate corr and print error if not zero 
               corr = flux_z(i,j,k,nspecies) + sumz 
               
               if(corr .gt. rhotot(i,j,k)*fraction_tolerance) then
                  !write(*,*) "Error: sum of z-fluxes = ", corr, " fluxes=", flux_z(i,j,k,:)
               end if
              
               ! correct z-flux for last species  
               flux_z(i,j,k,nspecies) = -sumz             
               total_corr = total_corr + abs(corr) 
 
            end do
         end do
      end do
      !write(*,*) "z flux correction = ", total_corr

    end subroutine correction_flux_3d

  end subroutine correction_flux

end module correction_flux_module
