module diffusive_m_fluxdiv_module

  use ml_layout_module
  use define_bc_module
  use bc_module
  use stag_applyop_module, only: stag_applyop_level

  implicit none

  private

  public :: diffusive_m_fluxdiv

contains

  subroutine diffusive_m_fluxdiv(mla,m_update,umac,eta,eta_ed,kappa,dx,the_bc_level)


    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(inout) :: m_update(:,:)
    type(multifab) , intent(in   ) ::     umac(:,:)
    type(multifab) , intent(in   ) ::      eta(:)
    type(multifab) , intent(in   ) ::   eta_ed(:,:)
    type(multifab) , intent(in   ) ::    kappa(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:)
    type(bc_level) , intent(in   ) :: the_bc_level(:)

    ! local variables
    integer :: nlevs,dm,i,n

    type(multifab) :: Lphi_fc(mla%nlevel,mla%dim)
    type(multifab) :: alpha_fc(mla%nlevel,mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"diffusive_m_fluxdiv")

    nlevs = mla%nlevel
    dm    = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(Lphi_fc(n,i),mla%la(n),1,0,i)
          call multifab_build_edge(alpha_fc(n,i),mla%la(n),1,0,i)
          ! set alpha to zero
          call setval(alpha_fc(n,i),0.d0,all=.true.)
       end do
    end do

    do n=1,nlevs

       ! compute -L(phi)
       ! we could compute +L(phi) but then we'd have to multiply beta and kappa by -1
       call stag_applyop_level(mla%la(n),the_bc_level(n),umac(n,:),Lphi_fc(n,:), &
                               alpha_fc(n,:),eta(n),eta_ed(n,:),kappa(n),dx(n,:))

       ! subtract -L(phi) to m_update
       do i=1,dm
          call multifab_sub_sub_c(m_update(n,i),1,Lphi_fc(n,i),1,1,0)
       end do

    end do

    do n=1,nlevs
       do i=1,dm
          call multifab_destroy(Lphi_fc(n,i))
          call multifab_destroy(alpha_fc(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine diffusive_m_fluxdiv

end module diffusive_m_fluxdiv_module
