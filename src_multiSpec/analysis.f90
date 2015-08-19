module analysis_module

  use bl_types
  use bl_constants_module
  use multifab_module
  use define_bc_module
  use ml_layout_module
  use init_lowmach_module
  use init_temp_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: total_volume
 
  implicit none

  private

  public :: print_errors,sum_mass,compute_cov

  contains

  subroutine print_errors(mla,rho,rho_exact,Temp,dx,time,the_bc_level)

    type(ml_layout) , intent(in)     :: mla
     type(multifab) , intent(in)     :: rho(:)            
     type(multifab) , intent(inout)  :: rho_exact(:)            
     type(multifab) , intent(inout)  :: Temp(:)            
     real(kind=dp_t), intent(in   )  :: dx(:,:)           
     real(kind=dp_t), intent(in   )  :: time 
     type(bc_level) , intent(in   )  :: the_bc_level(:)
 
     ! local variables
     real(kind=dp_t), dimension(nspecies) :: norm_inf,norm_l1,norm_l2
     real(kind=dp_t)                      :: norm_inf_tot,norm_l1_tot,norm_l2_tot
     integer                              :: i,n,nlevs,dm
     type(multifab)                       :: umac_tmp(mla%nlevel,mla%dim)

     type(bl_prof_timer), save :: bpt

     call build(bpt,"print_errors")

     nlevs = mla%nlevel
     dm = mla%dim

     do n=1,nlevs
        do i=1,dm
           call multifab_build_edge(umac_tmp(n,i),mla%la(n),1,1,i)
        end do
     end do

     ! calculate rho_exact
     call init_rho_and_umac(mla,rho_exact,umac_tmp,dx,time,the_bc_level) 
     call init_Temp(Temp,dx,time,the_bc_level) 
    
     ! substract the values 
     do n=1,nlevs
        call multifab_sub_sub_c(rho_exact(n),1,rho(n),1,nspecies,0)
     end do

     ! calculate norms for each species
     do n=1,nlevs
        do i=1, nspecies 
         
           ! Linf norm = max(diff_i)
           norm_inf(i) = multifab_norm_inf_c(rho_exact(n),i,1,all=.false.)

           ! L1 norm = 1/total_volume*sum(diff_i) 
           norm_l1(i) = multifab_sum_c(rho_exact(n),i,1,all=.false.)/total_volume

           ! L2 norm = sqrt{1/total_volume*sum(diff_i^2)} 
           norm_l2(i) = multifab_norm_l2_c(rho_exact(n),i,1,all=.false.)/sqrt(dble(total_volume))
     
        end do
     end do
     
     ! calculate total norms 
     norm_inf_tot = multifab_norm_inf_c(rho_exact(1),1,nspecies,all=.false.)
     norm_l1_tot = multifab_sum_c(rho_exact(1),1,nspecies,all=.false.)/total_volume
     norm_l2_tot = multifab_norm_l2_c(rho_exact(1),1,nspecies,all=.false.)/sqrt(dble(total_volume))
     
     ! print the norms
     if(.false.) then
        if (parallel_IOProcessor()) then 
            if(time.gt.2.99d0) then
            !if(time.gt.0.49d0) then
               write(*,*), time, norm_inf(1:nspecies), norm_l1(1:nspecies),norm_l2(1:nspecies),&
                           norm_inf_tot,norm_l1_tot,norm_l2_tot
            end if
        end if
     end if
     
     ! for checking analytic solution with Visit
     call init_rho_and_umac(mla,rho_exact,umac_tmp,dx,time,the_bc_level) 
     call init_Temp(Temp,dx,time,the_bc_level) 

     do n=1,nlevs
        do i=1,dm
           call multifab_destroy(umac_tmp(n,i))
        end do
     end do

     call destroy(bpt)

  end subroutine print_errors

  subroutine sum_mass(rho, step)

    type(multifab), intent(in   ) :: rho(:)
    integer,        intent(in   ) :: step

    ! local variables
    real(kind=dp_t) :: mass(nspecies)
    integer         :: i,n,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"sum_mass")

    nlevs = size(rho,1)
    
    do n=1,nlevs
       do i=1,nspecies
          mass(i) = multifab_sum_c(rho(n),i,1,all=.false.)
          if (parallel_IOProcessor()) then
             print*, step, ' sum of rho_i for i=',i,mass(i)
          end if
       end do
       if (parallel_IOProcessor()) then
          print*, step, ' sum of rho=',sum(mass(1:nspecies))
       end if
    end do

    call destroy(bpt)

  end subroutine sum_mass

  subroutine compute_cov(mla,rho,wit,wiwjt) 

    type(ml_layout), intent(in)     :: mla
    type(multifab),  intent(in)     :: rho(:)
    real(kind=dp_t), intent(inout)  :: wit(nspecies)
    real(kind=dp_t), intent(inout)  :: wiwjt(nspecies,nspecies)

    ! local variables
    integer :: n,i,j,dm,nlevs
    integer :: lo(rho(1)%dim), hi(rho(1)%dim)
    real(kind=dp_t), dimension(nspecies) :: cellW, cellW_procavg    
    real(kind=dp_t), dimension(nspecies,nspecies) :: cellWij, cellWij_procavg   
 
    ! pointer for rho 
    real(kind=dp_t), pointer  :: dp(:,:,:,:)  

    type(bl_prof_timer), save :: bpt

    call build(bpt, "compute_cov")

    dm     = mla%dim     ! dimensionality
    nlevs  = mla%nlevel 
 
    cellW_procavg   = 0.d0
    cellWij_procavg = 0.d0
    cellW           = 0.d0
    cellWij         = 0.d0
 
    ! loop over all boxes 
    do n=1,nlevs
       do i=1,nfabs(rho(n))
          dp => dataptr(rho(n),i)
          lo = lwb(get_box(rho(n),i))
          hi = upb(get_box(rho(n),i))

          select case(dm)
          case (2)
             call compute_cov_2d(dp(:,:,1,:),cellW,cellWij,lo,hi) 
          case (3)
             call compute_cov_3d(dp(:,:,:,:),cellW,cellWij,lo,hi) 
          end select
       end do
    end do
    
    ! average over all processors 
    call parallel_reduce(cellW_procavg, cellW, MPI_SUM)
    do i=1,nspecies
       do j=1, nspecies     
          call parallel_reduce(cellWij_procavg(i,j), cellWij(i,j), MPI_SUM)
       end do
    end do

    ! average over total_volume and calculate covW
    wit   = wit   + cellW_procavg/total_volume
    wiwjt = wiwjt + cellWij_procavg/total_volume 
 
    call destroy(bpt)

    end subroutine compute_cov
     
    subroutine compute_cov_2d(rho,cellW,cellWij,lo,hi)
 
       integer         :: lo(2), hi(2)
       real(kind=dp_t) :: rho(lo(1):,lo(2):,:) 
       real(kind=dp_t) :: cellW(nspecies)  
       real(kind=dp_t) :: cellWij(nspecies,nspecies)  
   
       ! local variables
       integer  :: i,j
              
       ! for specific box, now start loops over alloted cells   
       do j=lo(2), hi(2)
          do i=lo(1), hi(1)
             call compute_cov_local(rho(i,j,:),cellW,cellWij)
          end do
       end do

    end subroutine compute_cov_2d

    subroutine compute_cov_3d(rho,cellW,cellWij,lo,hi)
 
       integer          :: lo(3), hi(3)
       real(kind=dp_t)  :: rho(lo(1):,lo(2):,lo(3):,:) 
       real(kind=dp_t)  :: cellW(nspecies)  
       real(kind=dp_t)  :: cellWij(nspecies,nspecies)  
    
       ! local variables
       integer :: i,j,k
   
       ! for specific box, now start loops over alloted cells    
       do k=lo(3), hi(3)
          do j=lo(2), hi(2)
             do i=lo(1), hi(1)
                call compute_cov_local(rho(i,j,k,:),cellW,cellWij)
             end do
          end do
       end do

     end subroutine compute_cov_3d

     subroutine compute_cov_local(rho,cellW,cellWij)
 
       real(kind=dp_t)  :: rho(nspecies)  
       real(kind=dp_t)  :: cellW(nspecies)    
       real(kind=dp_t)  :: cellWij(nspecies,nspecies) 
    
       ! local variables
       real(kind=dp_t), dimension(nspecies)           :: W
       real(kind=dp_t), dimension(nspecies,nspecies)  :: Wij
       real(kind=dp_t)                                :: rhotot             
       integer                                        :: i,j

       type(bl_prof_timer), save :: bpt

       call build(bpt, "compute_cov_local")

       ! calculate total density inside each cell
       rhotot=0.d0 
       do i=1, nspecies  
          rhotot = rhotot + rho(i)
       end do         
  
       ! calculate mass fraction and sum over cell for each species
       do i=1, nspecies  
          W(i)     = rho(i)/rhotot
          cellW(i) = cellW(i) + W(i) ! compute this for average
       end do

       ! calculate Wij=wi*wj matrix and spatial sum over Wij
       do i=1,nspecies
          do j=1, i-1
             Wij(i,j) = W(i)*W(j)
             Wij(j,i) = Wij(i,j) 
             cellWij(i,j) = cellWij(i,j) + Wij(i,j)    
             cellWij(j,i) = cellWij(j,i) + Wij(j,i) 
          end do
          Wij(i,i) = W(i)*W(i)
          cellWij(i,i) = cellWij(i,i) + Wij(i,i) 
       end do 

       call destroy(bpt)
       
     end subroutine compute_cov_local 
 
end module analysis_module
