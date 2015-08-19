module mass_flux_utilities_module 

  use multifab_module
  use ml_layout_module
  use probin_common_module, only: k_B, molmass, rhobar
  use probin_multispecies_module, only: nspecies, use_lapack, fraction_tolerance, &
       is_ideal_mixture, inverse_type, is_nonisothermal, chi_iterations
  use matrix_utilities 
  use F95_LAPACK

  implicit none

  private

  public :: correct_rho_with_drho, &
            convert_cons_to_prim, &
            compute_rhotot, &
            compute_Gama, &
            compute_chi, &
            compute_rhoWchi, &
            compute_Lonsager, &
            compute_baro_coef

contains
  
  subroutine correct_rho_with_drho(mla,rho,drho)

   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(inout)  :: rho(:) 
   type(multifab) , intent(inout)  :: drho(:) 

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,ng_1,ng_2,dm,nlevs
 
   ! pointer for rho(nspecies), rhotot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for drho

   type(mfiter) :: mfi
   type(box) :: tilebox
   integer :: tlo(mla%dim), thi(mla%dim)

   type(bl_prof_timer), save :: bpt

   call build(bpt,"correct_rho_with_drho")

   dm    = mla%dim     ! dimensionality
   ng_1  = rho(1)%ng   ! number of ghost cells 
   ng_2  = drho(1)%ng
   nlevs = mla%nlevel  ! number of levels 
 
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi,dp1,dp2,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,rho(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(drho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))          
          select case(dm)
          case (2)
             call correct_rho_with_drho_2d(dp1(:,:,1,:),dp2(:,:,1,:),ng_1,ng_2,lo,hi,tlo,thi) 
          case (3)
             call correct_rho_with_drho_3d(dp1(:,:,:,:),dp2(:,:,:,:),ng_1,ng_2,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine correct_rho_with_drho

  subroutine correct_rho_with_drho_2d(rho,drho,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(2), ghi(2), ng_1, ng_2, tlo(2), thi(2)
    real(kind=dp_t)  ::  rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: drho(glo(1)-ng_2:,glo(2)-ng_2:,:) ! total density in each cell 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=tlo(2), thi(2)
       do i=tlo(1), thi(1)
         
         call correct_rho_with_drho_local(rho(i,j,:),drho(i,j,:))

       end do
    end do
 
  end subroutine correct_rho_with_drho_2d

  subroutine correct_rho_with_drho_3d(rho,drho,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_1, ng_2,tlo(:),thi(:)
    real(kind=dp_t)  ::  rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: drho(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! total density in each cell 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3), thi(3)
       do j=tlo(2), thi(2)
          do i=tlo(1), thi(1)

             call correct_rho_with_drho_local(rho(i,j,k,:),drho(i,j,k,:))

          end do
       end do
    end do
 
  end subroutine correct_rho_with_drho_3d

  subroutine correct_rho_with_drho_local(rho,drho)
 
    real(kind=dp_t), intent(inout) :: rho(nspecies)    ! density- last dim for #species
    real(kind=dp_t), intent(out)   :: drho(nspecies)   ! total density in each cell 
    
    ! local variables
    real(kind=dp_t)  :: rhotot_local
    real(kind=dp_t)  :: rho_old(nspecies)

    ! make a copy of the input rho's
    rho_old = rho
    
    ! total rho in the cell
    rhotot_local = sum(rho)

    ! make sure each rho is greater than +epsilon 
    ! to prevent division by zero in case species
    ! density, molar concentration or total density = 0.
    rho(1:nspecies) = max(rho(1:nspecies), fraction_tolerance*rhotot_local) 

    ! keep track of how much rho changed so we can reset it back later
    drho = rho - rho_old

  end subroutine correct_rho_with_drho_local 

  subroutine convert_cons_to_prim(mla,rho,rhotot,molarconc,molmtot)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in   )  :: rho(:) 
   type(multifab) , intent(inout)  :: rhotot(:) 
   type(multifab) , intent(inout)  :: molarconc(:) 
   type(multifab) , intent(inout)  :: molmtot(:) 

   ! local variables
   integer :: lo(rho(1)%dim), hi(rho(1)%dim)
   integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4
 
   ! pointer for rho(nspecies), rhotot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rhotot
   real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molar concentrations
   real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for total molar concentrations

   type(mfiter) :: mfi
   type(box) :: tilebox
   integer :: tlo(mla%dim), thi(mla%dim)

   type(bl_prof_timer), save :: bpt

   call build(bpt,"convert_cons_to_prim")

   dm    = mla%dim     ! dimensionality
   ng_1  = rho(1)%ng   ! number of ghost cells 
   ng_2  = rhotot(1)%ng
   ng_3  = molarconc(1)%ng
   ng_4  = molmtot(1)%ng
   nlevs = mla%nlevel  ! number of levels 

   if (ng_3 .ne. ng_4) then
      call bl_error('convert_cons_to_prim: ng for molarconc and molmass differ')
   end if
 
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi,dp1,dp2,dp3,dp4,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,molarconc(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

 !      do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(rhotot(n),i)
          dp3 => dataptr(molarconc(n),i)
          dp4 => dataptr(molmtot(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call compute_molconc_rhotot_2d(dp1(:,:,1,:),dp2(:,:,1,1), &
                                            dp3(:,:,1,:),dp4(:,:,1,1), &
                                            ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          case (3)
             call compute_molconc_rhotot_3d(dp1(:,:,:,:),dp2(:,:,:,1), &
                                            dp3(:,:,:,:),dp4(:,:,:,1), &
                                            ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine convert_cons_to_prim

  subroutine compute_molconc_rhotot_2d(rho,rhotot,molarconc,molmtot, &
                                       ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
 
    integer          :: glo(2), ghi(2), ng_1, ng_2, ng_3, ng_4
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::       rho(glo(1)-ng_1:,glo(2)-ng_1:,:)       ! density- last dim for #species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_2:,glo(2)-ng_2:)     ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_3:,glo(2)-ng_3:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_4:,glo(2)-ng_4:)     ! total molar mass 
        
    ! local variables
    integer          :: i,j
    
    ! for specific box, now start loops over alloted cells    
    do j=tlo(2), thi(2)
       do i=tlo(1), thi(1)
         
         call compute_molconc_rhotot_local(rho(i,j,:),rhotot(i,j), &
                                           molarconc(i,j,:),molmtot(i,j))

       end do
    end do
 
  end subroutine compute_molconc_rhotot_2d

  subroutine compute_molconc_rhotot_3d(rho,rhotot,molarconc,molmtot, &
                                       ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_1, ng_2, ng_3, ng_4
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::       rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:)   ! total molar mass 
    
    ! local variables
    integer          :: i,j,k
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3), thi(3)
       do j=tlo(2), thi(2)
          do i=tlo(1), thi(1)

             call compute_molconc_rhotot_local(rho(i,j,k,:),rhotot(i,j,k),&
                                               molarconc(i,j,k,:),molmtot(i,j,k))

          end do
       end do
    end do
 
  end subroutine compute_molconc_rhotot_3d

  subroutine compute_molconc_rhotot_local(rho,rhotot,molarconc,molmtot)
 
    real(kind=dp_t), intent(in)   :: rho(nspecies)       ! density- last dim for #species
    real(kind=dp_t), intent(out)  :: rhotot             ! total density in each cell 
    real(kind=dp_t), intent(out)  :: molarconc(nspecies) ! molar concentration
    real(kind=dp_t), intent(out)  :: molmtot             ! total molar mass 
    
    ! local variables
    integer          :: n
    real(kind=dp_t), dimension(nspecies) :: W            ! mass fraction w_i = rho_i/rho 
    real(kind=dp_t)  :: Sum_woverm, rhotot_local


    ! calculate total density inside each cell
    rhotot_local=0.d0 
    do n=1, nspecies  
       rhotot_local = rhotot_local + rho(n)
    end do         
    rhotot = rhotot_local

    ! calculate mass fraction and total molar mass (1/m=Sum(w_i/m_i))
    Sum_woverm=0.d0
    do n=1, nspecies  
       W(n) = rho(n)/rhotot
       Sum_woverm = Sum_woverm + W(n)/molmass(n)
    end do
    molmtot = 1.0d0/Sum_woverm 

    ! calculate molar concentrations in each cell (x_i=m*w_i/m_i) 
    do n=1, nspecies 
       molarconc(n) = molmtot*W(n)/molmass(n)
    end do
    

  end subroutine compute_molconc_rhotot_local 

  subroutine compute_rhotot(mla,rho,rhotot)
   
   type(ml_layout), intent(in   )  :: mla
   type(multifab) , intent(in   )  :: rho(:) 
   type(multifab) , intent(inout)  :: rhotot(:) 

   ! local variables
   integer :: lo(mla%dim), hi(mla%dim)
   integer :: n,i,dm,nlevs,ng_1,ng_2
 
   ! pointer for rho(nspecies), rhotot(1), molarconc(nspecies) 
   real(kind=dp_t), pointer        :: dp1(:,:,:,:)   ! for rho    
   real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rhotot

   type(mfiter) :: mfi
   type(box) :: tilebox
   integer :: tlo(mla%dim), thi(mla%dim)

   type(bl_prof_timer), save :: bpt

   call build(bpt, "compute_rhotot")

   dm    = mla%dim     ! dimensionality
   ng_1  = rho(1)%ng   ! number of ghost cells 
   ng_2  = rhotot(1)%ng
   nlevs = mla%nlevel  ! number of levels 
 
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi,dp1,dp2,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n),i)
          dp2 => dataptr(rhotot(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          
          select case(dm)
          case (2)
             call compute_rhotot_2d(dp1(:,:,1,:),dp2(:,:,1,1),&
                                    ng_1,ng_2,lo,hi,tlo,thi) 
          case (3)
             call compute_rhotot_3d(dp1(:,:,:,:),dp2(:,:,:,1),&
                                    ng_1,ng_2,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_rhotot

  subroutine compute_rhotot_2d(rho,rhotot,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(2), ghi(2), ng_1, ng_2, tlo(2), thi(2)
    real(kind=dp_t)  ::    rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: rhotot(glo(1)-ng_2:,glo(2)-ng_2:)   ! total density in each cell 
        
    ! local variables
    integer          :: i,j,n
    
    ! for specific box, now start loops over alloted cells    
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)

          ! calculate total density inside each cell
          rhotot(i,j)=0.d0
          do n=1, nspecies  
             rhotot(i,j) = rhotot(i,j) + rho(i,j,n)
          end do

       end do
    end do
 
  end subroutine compute_rhotot_2d

  subroutine compute_rhotot_3d(rho,rhotot,ng_1,ng_2,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_1, ng_2, tlo(3), thi(3)
    real(kind=dp_t)  ::    rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density- last dim for #species
    real(kind=dp_t)  :: rhotot(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:)   ! total density in each cell 
    
    ! local variables
    integer          :: i,j,k,n
    
    ! for specific box, now start loops over alloted cells    
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)

             ! calculate total density inside each cell
             rhotot(i,j,k)=0.d0
             do n=1,nspecies  
                rhotot(i,j,k) = rhotot(i,j,k) + rho(i,j,k,n)
             end do


          end do
       end do
    end do
 
  end subroutine compute_rhotot_3d

  subroutine compute_Gama(mla,rho,rhotot,molarconc,molmtot,Hessian,Gama)

    type(ml_layout), intent(in   )  :: mla
    type(multifab),  intent(in   )  :: rho(:) 
    type(multifab),  intent(in   )  :: rhotot(:) 
    type(multifab),  intent(in   )  :: molarconc(:) 
    type(multifab),  intent(in   )  :: molmtot(:) 
    type(multifab),  intent(in   )  :: Hessian(:)    ! Hessian matrix
    type(multifab),  intent(inout)  :: Gama(:)       ! Non-ideality coeficient 
 
    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,dm,nlevs,ng_0,ng_1,ng_2,ng_3,ng_4,ng_5

    ! assign pointers for multifabs to be passed
    real(kind=dp_t), pointer        :: dp0(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for Hessian 
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Gama 

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_Gama")

    dm    = mla%dim     ! dimensionality
    ng_0  = rho(1)%ng   ! number of ghost cells 
    ng_1  = rhotot(1)%ng
    ng_2  = molarconc(1)%ng
    ng_3  = molmtot(1)%ng
    ng_4  = Hessian(1)%ng
    ng_5  = Gama(1)%ng
    nlevs = mla%nlevel  ! number of levels 
  
   !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
   !$omp private(dp0,dp1,dp2,dp3,dp4,dp5,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, Gama(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,rho(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)
!       do i=1,nfabs(rho(n))
          dp0 => dataptr(rho(n),i)
          dp1 => dataptr(rhotot(n),i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(Hessian(n),i)
          dp5 => dataptr(Gama(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))          
          select case(dm)
          case (2)
             call compute_Gama_2d(dp0(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                                  dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:), &
                                  ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,lo,hi,tlo,thi) 
          case (3)
             call compute_Gama_3d(dp0(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                                  dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:), &
                                  ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_Gama
  
  subroutine compute_Gama_2d(rho,rhotot,molarconc,molmtot,Hessian,Gama, &
                             ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,glo,ghi,tlo,thi)

    integer          :: glo(2), ghi(2), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::       rho(glo(1)-ng_0:,glo(2)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_1:,glo(2)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_2:,glo(2)-ng_2:,:) ! molar concentration 
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_3:,glo(2)-ng_3:)   ! total molar mass 
    real(kind=dp_t)  ::   Hessian(glo(1)-ng_4:,glo(2)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(glo(1)-ng_5:,glo(2)-ng_5:,:) ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j

    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
       
          call compute_Gama_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                  molmtot(i,j),Hessian(i,j,:),Gama(i,j,:))

       end do
    end do
   
  end subroutine compute_Gama_2d

  subroutine compute_Gama_3d(rho,rhotot,molarconc,molmtot,Hessian,Gama, &
                             ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,glo,ghi,tlo,thi)
 
    integer          :: glo(3), ghi(3), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::       rho(glo(1)-ng_0:,glo(2)-ng_0:,glo(3)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::    rhotot(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! molar concentration; 
    real(kind=dp_t)  ::   molmtot(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:)   ! total molar mass 
    real(kind=dp_t)  ::   Hessian(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(glo(1)-ng_5:,glo(2)-ng_5:,glo(3)-ng_5:,:) ! last dimension for nspecies^2

    ! local varialbes
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)

             call compute_Gama_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                     molmtot(i,j,k),Hessian(i,j,k,:),Gama(i,j,k,:))
          end do
       end do
    end do
   
  end subroutine compute_Gama_3d

  subroutine compute_Gama_local(rho,rhotot,molarconc,molmtot,Hessian,Gama)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)        
    real(kind=dp_t), intent(in)   :: rhotot
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)
    real(kind=dp_t), intent(in)   :: molmtot
    real(kind=dp_t), intent(in)   :: Hessian(nspecies,nspecies)
    real(kind=dp_t), intent(out)  :: Gama(nspecies,nspecies)
 
    ! local variables
    integer                                       :: row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: I, X_xxT

    ! free the memory
    I=0; X_xxT=0;

    ! populate I and X_xxT, where X = molmtot*W*M^(-1) 
    do row=1, nspecies  
       do column=1, row-1
          I(row, column)    = 0.d0      
          I(column, row)    = I(row, column) ! symmetric
          
          if(.not. is_ideal_mixture) then
             X_xxT(row,column)   = -molarconc(row)*molarconc(column)  ! form x*transpose(x) off diagonals 
             X_xxT(column, row)  = X_xxT(row, column)                 ! symmetric
          end if
       end do
       
       I(row, row) = 1.d0        ! unit matrix for ideal mixture
       if(.not. is_ideal_mixture) then
          X_xxT(row,row) = molarconc(row) - molarconc(row)**2 
       end if
    end do
  
    ! compute Gama 
    Gama = I + matmul(X_xxT, Hessian)     
 
  end subroutine compute_Gama_local

  subroutine compute_chi(mla,rho,rhotot,molarconc,chi,D_bar,D_therm,Temp,zeta_by_Temp)
   
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: molarconc(:)
    type(multifab) , intent(inout)  :: chi(:) 
    type(multifab) , intent(in   )  :: D_bar(:) 
    type(multifab) , intent(in   )  :: D_therm(:) 
    type(multifab) , intent(in   )  :: Temp(:) 
    type(multifab) , intent(inout)  :: zeta_by_Temp(:) 

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,dm,nlevs,ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp0(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for D_bar
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Temp
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for zeta_by_Temp
    real(kind=dp_t), pointer        :: dp7(:,:,:,:)  ! for D_therm

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_chi")

    dm = mla%dim        ! dimensionality
    ng_0 = rho(1)%ng    ! number of ghost cells 
    ng_1 = rhotot(1)%ng
    ng_2 = molarconc(1)%ng
    ng_3 = chi(1)%ng
    ng_4 = D_bar(1)%ng
    ng_5 = Temp(1)%ng
    ng_6 = zeta_by_Temp(1)%ng
    ng_7 = D_therm(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(dp0,dp1,dp2,dp3,dp4,dp5,dp6,dp7,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,chi(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp0  => dataptr(rho(n), i)
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(molarconc(n), i)
          dp3 => dataptr(chi(n), i)
          dp4 => dataptr(D_bar(n), i)
          dp5 => dataptr(Temp(n), i)
          dp6 => dataptr(zeta_by_Temp(n), i)
          dp7 => dataptr(D_therm(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))          
          select case(dm)
          case (2)
             call compute_chi_2d(dp0(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),dp3(:,:,1,:),&
                                 dp4(:,:,1,:),dp5(:,:,1,1),dp6(:,:,1,:),dp7(:,:,1,:),&
                                 ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,lo,hi,tlo,thi) 
          case (3)
             call compute_chi_3d(dp0(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),dp3(:,:,:,:),&
                                 dp4(:,:,:,:),dp5(:,:,:,1),dp6(:,:,:,:),dp7(:,:,:,:),&
                                 ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_chi
 
  subroutine compute_chi_2d(rho,rhotot,molarconc,chi,D_bar,Temp,zeta_by_Temp,D_therm, &
                            ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,glo,ghi,tlo,thi)

    integer          :: glo(2), ghi(2), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7
    integer          :: tlo(2), thi(2)
    real(kind=dp_t)  ::          rho(glo(1)-ng_0:,glo(2)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::       rhotot(glo(1)-ng_1:,glo(2)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,:) ! molar concentration 
    real(kind=dp_t)  ::          chi(glo(1)-ng_3:,glo(2)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,:) ! MS diff-coefs 
    real(kind=dp_t)  ::         Temp(glo(1)-ng_5:,glo(2)-ng_5:)   ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(glo(1)-ng_6:,glo(2)-ng_6:,:) ! zeta/T
    real(kind=dp_t)  ::      D_therm(glo(1)-ng_7:,glo(2)-ng_7:,:) ! thermo diff-coefs 

    ! local variables
    integer          :: i,j,row
    real(kind=dp_t), dimension(nspecies,nspecies) :: chilocal

    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
    
          call compute_chi_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                 chi(i,j,:),D_bar(i,j,:),Temp(i,j),zeta_by_Temp(i,j,:),&
                                 D_therm(i,j,:))

          ! print chi for one cell 
          if(.false.) then 
          if(i.eq.32 .and. j.eq.16) then
            if(use_lapack) then 
              print*, 'print chi via inverse/p-inverse'
            else 
              print*, 'print chi via iterative methods'
            end if
            call set_Xij(chilocal, chi(i,j,:)) 
            do row=1, nspecies
               print*, chilocal(row, :)
            end do
          end if
          end if

       end do
    end do

  end subroutine compute_chi_2d

  subroutine compute_chi_3d(rho,rhotot,molarconc,chi,D_bar,Temp,zeta_by_Temp,D_therm, &
                            ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7,glo,ghi,tlo,thi)
   
    integer          :: glo(3), ghi(3), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,ng_7
    integer          :: tlo(3), thi(3)
    real(kind=dp_t)  ::          rho(glo(1)-ng_0:,glo(2)-ng_0:,glo(3)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::       rhotot(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  ::    molarconc(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:,:) ! molar concentration; 
    real(kind=dp_t)  ::          chi(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::        D_bar(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! SM diffusion constants 
    real(kind=dp_t)  ::         Temp(glo(1)-ng_5:,glo(2)-ng_5:,glo(3)-ng_5:)   ! Temperature 
    real(kind=dp_t)  :: zeta_by_Temp(glo(1)-ng_6:,glo(2)-ng_6:,glo(3)-ng_6:,:) ! zeta/T
    real(kind=dp_t)  ::      D_therm(glo(1)-ng_7:,glo(2)-ng_7:,glo(3)-ng_7:,:) ! thermo diffusion constants 
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_chi_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                    chi(i,j,k,:),D_bar(i,j,k,:),Temp(i,j,k),zeta_by_Temp(i,j,k,:),&
                                    D_therm(i,j,k,:))

          end do
       end do
    end do
   
  end subroutine compute_chi_3d

  subroutine compute_chi_local(rho,rhotot,molarconc,chi,D_bar,Temp,zeta_by_Temp,D_therm)
    
    real(kind=dp_t), intent(inout)  :: rho(nspecies)         
    real(kind=dp_t), intent(in)     :: rhotot               
    real(kind=dp_t), intent(inout)  :: molarconc(nspecies) 
    real(kind=dp_t), intent(out)    :: chi(nspecies,nspecies)   
    real(kind=dp_t), intent(in)     :: D_bar(nspecies,nspecies) 
    real(kind=dp_t), intent(in)     :: Temp
    real(kind=dp_t), intent(inout)  :: zeta_by_Temp(nspecies)
    real(kind=dp_t), intent(in)     :: D_therm(nspecies)

    ! local variables
    integer                         :: row,column
    real(kind=dp_t)                 :: Sum_knoti   

    ! vectors and matrices to be used by LAPACK 
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lambda
    real(kind=dp_t), dimension(nspecies)          :: W

    ! free up memory  
    Lambda        = 0.d0         
    W             = 0.d0

    ! compute Lambda_ij matrix and massfraction W_i = rho_i/rho; molarconc is 
    ! expressed in terms of molmtot,mi,rhotot etc. 
    do row=1, nspecies  
       do column=1, row-1
          Lambda(row, column) = -molarconc(row)*molarconc(column)/D_bar(row,column)
          Lambda(column, row) = Lambda(row, column) 
       end do
       W(row) = rho(row)/rhotot
    end do

    ! compute Lambda_ii
    do row=1, nspecies
       Sum_knoti = 0.d0
       do column=1, nspecies
          if(column.ne.row) then
             Sum_knoti = Sum_knoti - Lambda(row,column)
          end if
          Lambda(row,row) = Sum_knoti
       end do
    end do

    ! compute zeta_by_Temp for thermodiffusion
    if(is_nonisothermal) then
       do row=1, nspecies
          Sum_knoti = 0.d0
          do column=1, nspecies
             if(column.ne.row) then
                Sum_knoti = Sum_knoti + Lambda(row,column)*(D_therm(row)-D_therm(column))
             end if
             zeta_by_Temp(row) = Sum_knoti/Temp
          end do
       end do
    end if
    
    ! compute chi either selecting inverse/pseudoinverse or iterative methods 
    if(use_lapack) then
       call compute_chi_lapack(Lambda(:,:),chi(:,:),W(:))
    else
       call Dbar2chi_iterative(chi_iterations,D_bar(:,:),molarconc(:),chi(:,:)) 
    end if

  end subroutine compute_chi_local

  subroutine compute_chi_lapack(Lambda,chi,W)
    
    real(kind=dp_t)  :: Lambda(nspecies,nspecies)
    real(kind=dp_t)  :: chi(nspecies,nspecies)
    real(kind=dp_t)  :: W(nspecies)
 
    ! local variables
    integer          :: row,column,info
    real(kind=dp_t)  :: alpha    

    ! vectors and matrices to be used by LAPACK
    real(kind=dp_t), dimension(nspecies,nspecies) :: Sdag,chilocal
    real(kind=dp_t), dimension(nspecies,nspecies) :: U, UT, V, VT
    real(kind=dp_t), dimension(nspecies)          :: S, work 
    integer,         dimension(nspecies)          :: ipiv

    ! free up the memory  
    Sdag     = 0.d0
    U        = 0.d0
    UT       = 0.d0
    V        = 0.d0
    VT       = 0.d0
    S        = 0.d0
    work     = 0.d0
    alpha    = 0.d0
    chilocal = 0.d0
 
    ! calculate trace(Lambda)
    alpha = 0.d0
    do row=1, nspecies
       alpha = alpha + Lambda(row,row)
    end do
 
    ! calculate Lambda + alpha*W*WT (Equation 6) 
    do row=1, nspecies
       do column=1, nspecies
          chilocal(row,column) = alpha*W(row)*W(column) + Lambda(row,column)
       end do
    end do

    !===============================================================          
    ! select LAPACK inversion type, 1=inverse, 2=pseudo inverse 
    !===============================================================          
    select case(inverse_type) 
           
    case(1)
    !==========================================================
    ! Using Inverse 
    !==========================================================
 
    ! compute chilocal inverse
    call dgetrf(nspecies, nspecies, chilocal, nspecies, ipiv, info) 
    call dgetri(nspecies, chilocal, nspecies, ipiv, work, nspecies, info) 

    ! populate chi with B^(-1)
    chi = chilocal   
 
    case(2) 
    !==========================================================
    ! Using pseudoinverse 
    !==========================================================

    ! SVD decomposition of chilocal = U * S * VTranspose; note that chilocal 
    ! is changed. also V=(VT)T, UT = (U)T are needed for pseudoinverse of chilocal.
    call la_gesvd(chilocal, S, U, VT)
    V = transpose(VT)
    UT = transpose(U)
   
    ! populate diagonal matrix Sdag = 1/S with diagonal=0 below a chosen tolerance
    do row=1, nspecies
       do column=1,nspecies
          Sdag(row,column) = 0.0d0
       end do
       
       if(S(row).gt.fraction_tolerance*sum(S)) then 
          Sdag(row,row) = 1.0d0/S(row)
       else
          Sdag(row,row) = 0.0d0
       end if 
    end do

    ! compute chi = V*Sdag*UT, the pseudoinverse of chilocal 
    chi = matmul(V, matmul(Sdag, UT))

    end select
    !===============================================================          
 
    ! compute chi as equation (6) 
    do row=1, nspecies
       do column=1, nspecies
          chi(row,column) = chi(row,column) - 1.0d0/alpha
       end do
    end do
          
  end subroutine compute_chi_lapack

  subroutine compute_Lonsager(mla,rho,rhotot,molarconc,molmtot,chi,Gama,Lonsager)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: molarconc(:) 
    type(multifab) , intent(in   )  :: molmtot(:) 
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(in   )  :: Gama(:) 
    type(multifab) , intent(inout)  :: Lonsager(:)

    ! local variables
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    integer :: n,i,dm,nlevs,ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp0(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for molarconc
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for molmtot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp5(:,:,:,:)  ! for Gama
    real(kind=dp_t), pointer        :: dp6(:,:,:,:)  ! for Lonsager 

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_Lonsager")

    dm = mla%dim        ! dimensionality
    ng_0 = rho(1)%ng    ! number of ghost cells 
    ng_1 = rhotot(1)%ng
    ng_2 = molarconc(1)%ng
    ng_3 = molmtot(1)%ng
    ng_4 = chi(1)%ng
    ng_5 = Gama(1)%ng
    ng_6 = Lonsager(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp0 => dataptr(rho(n), i)
          dp1 => dataptr(rhotot(n), i)
          dp2 => dataptr(molarconc(n),i)
          dp3 => dataptr(molmtot(n),i)
          dp4 => dataptr(chi(n), i)
          dp5 => dataptr(Gama(n), i)
          dp6 => dataptr(Lonsager(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_Lonsager_2d(dp0(:,:,1,:),dp1(:,:,1,1),dp2(:,:,1,:),&
                                      dp3(:,:,1,1),dp4(:,:,1,:),dp5(:,:,1,:),dp6(:,:,1,:),&
                                      ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,lo,hi) 
          case (3)
             call compute_Lonsager_3d(dp0(:,:,:,:),dp1(:,:,:,1),dp2(:,:,:,:),&
                                      dp3(:,:,:,1),dp4(:,:,:,:),dp5(:,:,:,:),dp6(:,:,:,:),&
                                      ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,lo,hi) 
          end select
       end do
    end do
 
    call destroy(bpt)

  end subroutine compute_Lonsager
  
  subroutine compute_Lonsager_2d(rho,rhotot,molarconc,molmtot,chi,Gama,Lonsager, &
                                 ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,lo,hi)
  
    integer          :: lo(2), hi(2), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6
    real(kind=dp_t)  ::       rho(lo(1)-ng_0:,lo(2)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::    rhotot(lo(1)-ng_1:,lo(2)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng_2:,lo(2)-ng_2:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(lo(1)-ng_3:,lo(2)-ng_3:)   ! total molar mass 
    real(kind=dp_t)  ::       chi(lo(1)-ng_4:,lo(2)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(lo(1)-ng_5:,lo(2)-ng_5:,:) ! non-ideality coeficient 
    real(kind=dp_t)  ::  Lonsager(lo(1)-ng_6:,lo(2)-ng_6:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row
    real(kind=dp_t), dimension(nspecies,nspecies) :: Lonsager_local 
 
    ! for specific box, now start loops over alloted cells 
    do j=lo(2)-ng_6,hi(2)+ng_6
       do i=lo(1)-ng_6,hi(1)+ng_6
        
          call compute_Lonsager_local(rho(i,j,:),rhotot(i,j),molarconc(i,j,:),&
                                      molmtot(i,j),chi(i,j,:),Gama(i,j,:),Lonsager(i,j,:))

          if(.false.) then
          if(i.eq.7 .and. j.eq.14) then
          call set_Xij(Lonsager_local, Lonsager(i,j,:)) 
            write(*,*), "Lonsager"
            do row=1, nspecies
               write(*,*), Lonsager_local(row,:) 
               !write(*,*), sum(Lonsager_local(row,:))
            end do
            do row=1, nspecies
               !write(*,*), sum(Lonsager_local(:,row))
            end do
          end if
          end if

       end do
    end do

  end subroutine compute_Lonsager_2d

  subroutine compute_Lonsager_3d(rho,rhotot,molarconc,molmtot,chi,Gama,Lonsager, &
                                 ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6,lo,hi)

    integer          :: lo(3), hi(3), ng_0,ng_1,ng_2,ng_3,ng_4,ng_5,ng_6
    real(kind=dp_t)  ::       rho(lo(1)-ng_0:,lo(2)-ng_0:,lo(3)-ng_0:,:) ! density; last dimension for species
    real(kind=dp_t)  ::    rhotot(lo(1)-ng_1:,lo(2)-ng_1:,lo(3)-ng_1:)   ! total density in each cell 
    real(kind=dp_t)  :: molarconc(lo(1)-ng_2:,lo(2)-ng_2:,lo(3)-ng_2:,:) ! molar concentration
    real(kind=dp_t)  ::   molmtot(lo(1)-ng_3:,lo(2)-ng_3:,lo(3)-ng_3:)   ! total molar mass 
    real(kind=dp_t)  ::       chi(lo(1)-ng_4:,lo(2)-ng_4:,lo(3)-ng_4:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  ::      Gama(lo(1)-ng_5:,lo(2)-ng_5:,lo(3)-ng_5:,:) ! non-ideality coeficient 
    real(kind=dp_t)  ::  Lonsager(lo(1)-ng_6:,lo(2)-ng_6:,lo(3)-ng_6:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=lo(3)-ng_6,hi(3)+ng_6
       do j=lo(2)-ng_6,hi(2)+ng_6
          do i=lo(1)-ng_6,hi(1)+ng_6
       
             call compute_Lonsager_local(rho(i,j,k,:),rhotot(i,j,k),molarconc(i,j,k,:),&
                                         molmtot(i,j,k),chi(i,j,k,:),Gama(i,j,k,:),Lonsager(i,j,k,:))
              
         end do
      end do
    end do
   
  end subroutine compute_Lonsager_3d

subroutine compute_Lonsager_local(rho,rhotot,molarconc,molmtot,chi,Gama,Lonsager)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: rhotot                  
    real(kind=dp_t), intent(in)   :: molarconc(nspecies)      ! molar concentration
    real(kind=dp_t), intent(in)   :: molmtot                  ! total molar mass 
    real(kind=dp_t), intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
    real(kind=dp_t), intent(in)   :: Gama(nspecies,nspecies)        
    real(kind=dp_t), intent(out)  :: Lonsager(nspecies,nspecies) 

    ! local variables
    integer                              :: row,column,info
    real(kind=dp_t), dimension(nspecies) :: W 
    real(kind=dp_t)                      :: rcond 
  
    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_Lonsager_local")

    ! compute massfraction W_i = rho_i/rho; 
    do row=1, nspecies  
       W(row) = rho(row)/rhotot
    end do

    ! compute Onsager matrix L
    do column=1, nspecies
       do row=1, nspecies
          Lonsager(row, column) = molmtot*rhotot*W(row)*chi(row,column)*W(column)/k_B
       end do
    end do

    ! compute cell-centered Cholesky factor of Lonsager
    if(use_lapack) then
       
       call dpotrf_f95(Lonsager,'L', rcond, 'I', info)
    
       ! remove all upper-triangular entries and NXN entry that lapack doesn't set to zero 
       do row=1, nspecies
          do column=row+1, nspecies
             Lonsager(row, column) = 0.0d0          
          end do
       end do    
       Lonsager(nspecies, nspecies) = 0.0d0          
    
    else
       call choldc(Lonsager,nspecies)   
    end if

    call destroy(bpt)

  end subroutine compute_Lonsager_local

  subroutine compute_rhoWchi(mla,rho,rhotot,chi,rhoWchi)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: chi(:) 
    type(multifab) , intent(inout)  :: rhoWchi(:) 

    ! local variables
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4
 
    ! pointer for rho(nspecies), molarconc(nspecies) 
    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for rho    
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for chi
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for rhoWchi

    type(mfiter) :: mfi
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_rhoWchi")

    dm = mla%dim        ! dimensionality
    ng_1 = rho(1)%ng    ! number of ghost cells 
    ng_2 = rhotot(1)%ng
    ng_3 = chi(1)%ng
    ng_4 = rhoWchi(1)%ng
    nlevs = mla%nlevel  ! number of levels 
 
    !$omp parallel private(n,i,mfi,tilebox,tlo,thi) &
    !$omp private(dp1,dp2,dp3,dp4,lo,hi)

    ! loop over all boxes 
    do n=1,nlevs
       call mfiter_build(mfi, rho(n), tiling=.true.)

       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_growntilebox(mfi,rhoWchi(n)%ng)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

!       do i=1,nfabs(rho(n))
          dp1 => dataptr(rho(n), i)
          dp2 => dataptr(rhotot(n), i)
          dp3 => dataptr(chi(n), i)
          dp4 => dataptr(rhoWchi(n), i)
          lo  =  lwb(get_box(rho(n), i))
          hi  =  upb(get_box(rho(n), i))
          
          select case(dm)
          case (2)
             call compute_rhoWchi_2d(dp1(:,:,1,:),dp2(:,:,1,1),dp3(:,:,1,:),dp4(:,:,1,:), &
                                     ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          case (3)
             call compute_rhoWchi_3d(dp1(:,:,:,:),dp2(:,:,:,1),dp3(:,:,:,:),dp4(:,:,:,:), &
                                     ng_1,ng_2,ng_3,ng_4,lo,hi,tlo,thi) 
          end select
       end do
    end do
    !$omp end parallel

    call destroy(bpt)

  end subroutine compute_rhoWchi
  
  subroutine compute_rhoWchi_2d(rho,rhotot,chi,rhoWchi,ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)
  
    integer          :: glo(2), ghi(2), ng_1,ng_2,ng_3,ng_4,tlo(2),thi(2)
    real(kind=dp_t)  ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,:) ! density; last dimension for species
    real(kind=dp_t)  ::  rhotot(glo(1)-ng_2:,glo(2)-ng_2:)   ! total density in each cell
    real(kind=dp_t)  ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,:) ! last dimension for nspecies^2

    ! local variables
    integer          :: i,j,row,column
    real(kind=dp_t), dimension(nspecies,nspecies) :: rhoWchiloc 
  
    ! for specific box, now start loops over alloted cells 
    do j=tlo(2),thi(2)
       do i=tlo(1),thi(1)
        
          call compute_rhoWchi_local(rho(i,j,:),rhotot(i,j),chi(i,j,:),rhoWchi(i,j,:))

          if(.false.) then
          if(i.eq.7 .and. j.eq.14) then
          call set_Xij(rhoWchiloc, rhoWchi(i,j,:))
             do row=1, nspecies
                do column=1, nspecies
                   print*, rhoWchiloc(row,column) 
                end do
                print*, '' 
             end do
          end if
          end if

       end do
    end do

  end subroutine compute_rhoWchi_2d

  subroutine compute_rhoWchi_3d(rho,rhotot,chi,rhoWchi,ng_1,ng_2,ng_3,ng_4,glo,ghi,tlo,thi)

    integer          :: glo(3), ghi(3), ng_1,ng_2,ng_3,ng_4,tlo(3),thi(3)
    real(kind=dp_t)  ::     rho(glo(1)-ng_1:,glo(2)-ng_1:,glo(3)-ng_1:,:) ! density; last dimension for species
    real(kind=dp_t)  ::  rhotot(glo(1)-ng_2:,glo(2)-ng_2:,glo(3)-ng_2:)   ! total density in each cell 
    real(kind=dp_t)  ::     chi(glo(1)-ng_3:,glo(2)-ng_3:,glo(3)-ng_3:,:) ! last dimension for nspecies^2
    real(kind=dp_t)  :: rhoWchi(glo(1)-ng_4:,glo(2)-ng_4:,glo(3)-ng_4:,:) ! last dimension for nspecies^2
    
    ! local variables
    integer          :: i,j,k

    ! for specific box, now start loops over alloted cells 
    do k=tlo(3),thi(3)
       do j=tlo(2),thi(2)
          do i=tlo(1),thi(1)
       
             call compute_rhoWchi_local(rho(i,j,k,:),rhotot(i,j,k),chi(i,j,k,:),rhoWchi(i,j,k,:))
              
         end do
      end do
    end do
   
  end subroutine compute_rhoWchi_3d
  
  subroutine compute_rhoWchi_local(rho,rhotot,chi,rhoWchi)
   
    real(kind=dp_t), intent(in)   :: rho(nspecies)            
    real(kind=dp_t), intent(in)   :: rhotot                  
    real(kind=dp_t), intent(in)   :: chi(nspecies,nspecies)   ! rank conversion done 
    real(kind=dp_t), intent(out)  :: rhoWchi(nspecies,nspecies) 
 
    ! local variables
    integer                              :: row,column
    real(kind=dp_t), dimension(nspecies) :: W !,chiw 

    ! compute massfraction W_i = rho_i/rho; 
    do row=1, nspecies  
       W(row) = rho(row)/rhotot
    end do

    ! populate -rho*W*chi = -rho_i*chi
    do row=1, nspecies
       do column=1, nspecies
          rhoWchi(row,column) = -rho(row)*chi(row,column)  
       end do
    end do

  end subroutine compute_rhoWchi_local

  subroutine set_Xij(Xout_ij, Xin_ij)
        
    real(kind=dp_t), dimension(nspecies,nspecies), intent(in)  :: Xin_ij
    real(kind=dp_t), dimension(nspecies,nspecies), intent(out) :: Xout_ij  
   
    type(bl_prof_timer), save :: bpt

    call build(bpt, "set_Xij")

    ! reshape array into matrix without doing index algebra
    Xout_ij = Xin_ij
  
    call destroy(bpt)

  end subroutine set_Xij

  subroutine compute_baro_coef(mla,baro_coef,rho,rhotot,Temp)
 
    type(ml_layout), intent(in   )  :: mla
    type(multifab) , intent(inout)  :: baro_coef(:)
    type(multifab) , intent(in   )  :: rho(:)
    type(multifab) , intent(in   )  :: rhotot(:) 
    type(multifab) , intent(in   )  :: Temp(:) 

    ! local
    integer :: lo(mla%dim), hi(mla%dim)
    integer :: n,i,dm,nlevs,ng_1,ng_2,ng_3,ng_4

    real(kind=dp_t), pointer        :: dp1(:,:,:,:)  ! for baro_coef
    real(kind=dp_t), pointer        :: dp2(:,:,:,:)  ! for rho
    real(kind=dp_t), pointer        :: dp3(:,:,:,:)  ! for rhotot
    real(kind=dp_t), pointer        :: dp4(:,:,:,:)  ! for Temp

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_baro_coef")

    dm = mla%dim
    nlevs = mla%nlevel

    ng_1 = baro_coef(1)%ng
    ng_2 = rho(1)%ng
    ng_3 = rhotot(1)%ng
    ng_4 = Temp(1)%ng

    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp1 => dataptr(baro_coef(n),i)
          dp2 => dataptr(rho(n),i)
          dp3 => dataptr(rhotot(n),i)
          dp4 => dataptr(Temp(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))
          select case(dm)
          case (2)
             call compute_baro_coef_2d(dp1(:,:,1,:),dp2(:,:,1,:),dp3(:,:,1,1),dp4(:,:,1,1), &
                                        ng_1,ng_2,ng_3,ng_4,lo,hi) 
          case (3)
             call bl_error("compute_baro_coef_3d not written yet")
          end select
       end do
    end do

    call destroy(bpt)

  end subroutine compute_baro_coef

  subroutine compute_baro_coef_2d(baro_coef,rho,rhotot,Temp,ng_1,ng_2,ng_3,ng_4,lo,hi)
 
    integer        , intent(in   ) :: lo(2), hi(2), ng_1, ng_2, ng_3, ng_4
    real(kind=dp_t), intent(inout) :: baro_coef(lo(1)-ng_1:,lo(2)-ng_1:,:)
    real(kind=dp_t), intent(in   ) ::       rho(lo(1)-ng_2:,lo(2)-ng_2:,:)
    real(kind=dp_t), intent(in   ) ::    rhotot(lo(1)-ng_3:,lo(2)-ng_3:)
    real(kind=dp_t), intent(in   ) ::      Temp(lo(1)-ng_4:,lo(2)-ng_4:)
    
    ! local variables
    integer :: i,j,comp
    real(kind=dp_t) :: n

    do j=lo(2)-ng_1,hi(2)+ng_1
       do i=lo(1)-ng_1,hi(1)+ng_1

          n = 0.d0
          do comp=1,nspecies
             n = n + rho(i,j,comp)/molmass(comp)
          end do

          do comp=1,nspecies
             baro_coef(i,j,comp) = rho(i,j,comp)/rhobar(comp) / (n*k_B*Temp(i,j))
          end do

       end do
    end do

  end subroutine compute_baro_coef_2d
  
end module mass_flux_utilities_module
