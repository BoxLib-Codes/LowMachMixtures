module convert_rhoc_to_c_module

  use multifab_module
  use ml_layout_module
  use define_bc_module
  use bc_module
  use multifab_physbc_module
  use probin_common_module, only: rhobar
  use probin_multispecies_module, only: nspecies

  implicit none

  private

  public :: convert_rhoc_to_c, fill_c_ghost_cells
  
contains

  subroutine convert_rhoc_to_c(mla,rho,rhotot,conc,rho_to_c)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: rho(:)
    type(multifab) , intent(in   ) :: rhotot(:)
    type(multifab) , intent(inout) :: conc(:)
    logical        , intent(in   ) :: rho_to_c

    ! local
    integer :: n,nlevs,i

    type(bl_prof_timer), save :: bpt

    call build(bpt, "convert_rhoc_to_c")

    nlevs = mla%nlevel

    if (rho_to_c) then

       ! rho to conc - NO GHOST CELLS
       do n=1,nlevs
          call multifab_copy_c(conc(n),1,rho(n),1,nspecies,0)
          do i=1,nspecies
             call multifab_div_div_c(conc(n),i,rhotot(n),1,1,0)
          end do
       end do

    else

       ! conc to rho - VALID + GHOST (CAN CHANGE TO DO ONLY GHOST TO SAVE COMPUTATION)
       do n=1,nlevs
          call multifab_copy_c(rho(n),1,conc(n),1,nspecies,rho(n)%ng)
          do i=1,nspecies
             call multifab_mult_mult_c(rho(n),i,rhotot(n),1,1,rho(n)%ng)
          end do
       end do

    end if

    call destroy(bpt)

  end subroutine convert_rhoc_to_c

  subroutine fill_c_ghost_cells(mla,conc,dx,the_bc_tower)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: conc(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_tower) , intent(in   ) :: the_bc_tower

    ! local
    integer :: n,nlevs

    type(bl_prof_timer), save :: bpt

    call build(bpt,"fill_c_ghost_cells")

    nlevs = mla%nlevel

    do n=1,nlevs
       ! fill ghost cells for two adjacent grids including periodic boundary ghost cells
       call multifab_fill_boundary(conc(n))
       ! fill non-periodic domain boundary ghost cells
       call multifab_physbc(conc(n),1,c_bc_comp,nspecies,the_bc_tower%bc_tower_array(n), &
                            dx_in=dx(n,:))
    end do

    call destroy(bpt)

  end subroutine fill_c_ghost_cells

end module convert_rhoc_to_c_module
