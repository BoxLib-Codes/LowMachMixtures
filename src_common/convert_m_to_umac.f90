module convert_m_to_umac_module

  use multifab_module
  use ml_layout_module

  implicit none

  private

  public :: convert_m_to_umac
  
contains

  subroutine convert_m_to_umac(mla,s_fc,m,umac,m_to_umac)
    
    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) ::    s_fc(:,:)
    type(multifab) , intent(inout) ::    m(:,:)  
    type(multifab) , intent(inout) :: umac(:,:)
    logical        , intent(in   ) :: m_to_umac

    ! local
    integer :: n,i,dm,nlevs

    type(bl_prof_timer),save :: bpt

    call build(bpt,"convert_m_to_umac")

    dm = mla%dim
    nlevs = mla%nlevel

    if (m_to_umac) then

       ! compute umac = m / rho - NO GHOST CELLS
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(umac(n,i), 1, m(n,i), 1, 1, 0)
             call multifab_div_div_c(umac(n,i), 1, s_fc(n,i), 1, 1, 0)
          end do
       end do

    else

       ! compute m = rho * umac - INCLUDING GHOST CELLS
       do n=1,nlevs
          do i=1,dm
             call multifab_copy_c(m(n,i), 1, umac(n,i), 1, 1, m(n,i)%ng)
             call multifab_mult_mult_c(m(n,i), 1, s_fc(n,i), 1, 1, m(n,i)%ng)
          end do
       end do

    end if

    call destroy(bpt)

  end subroutine convert_m_to_umac

end module convert_m_to_umac_module
