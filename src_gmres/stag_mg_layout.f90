module stag_mg_layout_module

  use ml_layout_module
  use probin_common_module, only: n_cells
  use probin_gmres_module, only: stag_mg_max_bottom_nlevels, stag_mg_minwidth

  implicit none

  private

  public :: stag_mg_layout_build, stag_mg_layout_destroy, compute_nlevs_mg, &
       la_mg_normal, la_mg_fancy, mla_fancy
  
  type(layout), save, allocatable :: la_mg_normal(:)
  type(layout), save, allocatable :: la_mg_fancy(:)

  type(ml_layout), save :: mla_fancy

contains

  subroutine stag_mg_layout_build(mla)

    type(ml_layout), intent(in   ) :: mla

    ! hold the problem domain and boxarray at level 1 as well as the current multigrid level
    type(box)      :: pd_base,pd
    type(boxarray) :: ba_base,ba

    integer :: nlevs_mg,n,dm

    type(ml_boxarray) :: mba_fancy

    integer :: lo(mla%dim), hi(mla%dim), bottom_box_size
    type(box) :: bx

    dm = mla%dim
  
    ! get the problem domain from level 1
    pd_base = ml_layout_get_pd(mla,1)

    ! get boxarray from level 1
    ba_base = get_boxarray(mla%la(1))

    ! compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    ! smallest dimension of the smallest grid at the coarsest multigrid level
    call compute_nlevs_mg(nlevs_mg,ba_base)

    allocate(la_mg_normal(nlevs_mg))

    do n=1,nlevs_mg

       ! create the problem domain for this multigrid level
       pd = coarsen(pd_base,2**(n-1))

       ! create the boxarray for this multigrid level
       call boxarray_build_copy(ba,ba_base)
       call boxarray_coarsen(ba,2**(n-1))

       ! sanity check to make sure level 1 boxarrays match
       if (n .eq. 1 .and. (.not. boxarray_same_q(mla%mba%bas(1),ba) ) ) then
          call print(ba)
          call print(mla%mba%bas(1))
          call bl_error("Finest multigrid level boxarray and coarsest problem boxarrays do not match")
       end if

       ! build the layout, la
       ! force the same processor assignments as mla%la(1).  We can do this since there
       ! are the same number of boxes in the same order in physical space
       if (n .eq. 1) then
          la_mg_normal(1) = mla%la(1)
       else
          call layout_build_ba(la_mg_normal(n),ba,pd,mla%pmask,explicit_mapping=get_proc(mla%la(1)))
       end if

       ! don't need this anymore - free up memory
       call destroy(ba)

    end do

    !!!!!!!!!!!!!!!!!!!!!!!
    ! fancy bottom solver

    ! tell mba how many levels and dmensionality of problem
    call ml_boxarray_build_n(mba_fancy,1,mla%dim)
          
    ! create a box containing number of cells of current bottom solve
    lo(1:dm) = 0
    hi(1:dm) = n_cells(1:dm) / 2**(nlevs_mg-1) - 1
    
    bx = make_box(lo,hi)

    ! tell mba about the problem domain
    mba_fancy%pd(1) = bx

    ! initialize the boxarray at level 1 to be one single box
    call boxarray_build_bx(mba_fancy%bas(1),bx)

    ! chop up the box to respect stag_mg_max_bottom_nlevels
    bottom_box_size = maxval(hi(1:dm))+1
    bottom_box_size = min(bottom_box_size,stag_mg_minwidth*2**(stag_mg_max_bottom_nlevels-1))
    call boxarray_maxsize(mba_fancy%bas(1),bottom_box_size)

    ! build the ml_layout, mla_fancy
    call ml_layout_build(mla_fancy,mba_fancy,mla%pmask)

    if (parallel_IOProcessor()) then
       print*,'Building layouts for Fancy Bottom Solver'
       print*,'# of grids of old bottom solve:',nboxes(mla%la(1))
       print*,'# of grids of new bottom solve:',nboxes(mla_fancy%la(1))
    end if

    ! don't need this anymore - free up memory
    call destroy(mba_fancy)
  
    ! get the problem domain from level 1
    pd_base = ml_layout_get_pd(mla_fancy,1)

    ! get boxarray from level 1
    ba_base = get_boxarray(mla_fancy%la(1))

    ! compute the number of multigrid levels assuming stag_mg_minwidth is the length of the
    ! smallest dimension of the smallest grid at the coarsest multigrid level
    call compute_nlevs_mg(nlevs_mg,ba_base)

    allocate(la_mg_fancy(nlevs_mg))

    do n=1,nlevs_mg

       ! create the problem domain for this multigrid level
       pd = coarsen(pd_base,2**(n-1))

       ! create the boxarray for this multigrid level
       call boxarray_build_copy(ba,ba_base)
       call boxarray_coarsen(ba,2**(n-1))

       ! sanity check to make sure level 1 boxarrays match
       if (n .eq. 1 .and. (.not. boxarray_same_q(mla_fancy%mba%bas(1),ba) ) ) then
          call print(ba)
          call print(mla_fancy%mba%bas(1))
          call bl_error("Finest multigrid level boxarray and coarsest problem boxarrays do not match")
       end if

       ! build the layout, la
       ! force the same processor assignments as mla_fancy%la(1).  We can do this since there
       ! are the same number of boxes in the same order in physical space
       if (n .eq. 1) then
          la_mg_fancy(1) = mla_fancy%la(1)
       else
          call layout_build_ba(la_mg_fancy(n),ba,pd,mla_fancy%pmask,explicit_mapping=get_proc(mla_fancy%la(1)))
       end if

       ! don't need this anymore - free up memory
       call destroy(ba)

    end do

  end subroutine stag_mg_layout_build

  subroutine stag_mg_layout_destroy()

    integer :: n

    do n=2,size(la_mg_normal)
       call destroy(la_mg_normal(n))
    end do

    call destroy(mla_fancy)
    do n=2,size(la_mg_fancy)
       call destroy(la_mg_fancy(n))
    end do

    deallocate(la_mg_normal)
    deallocate(la_mg_fancy)

  end subroutine stag_mg_layout_destroy

  ! compute the number of multigrid levels assuming minwidth is the length of the
  ! smallest dimension of the smallest grid at the coarsest multigrid level
  subroutine compute_nlevs_mg(nlevs_mg,ba)

    integer       , intent(inout) :: nlevs_mg
    type(boxarray), intent(in   ) :: ba

    ! local
    integer :: i,d,dm,rdir,temp
    integer :: lo(get_dim(ba)),hi(get_dim(ba)),length(get_dim(ba))

    type(bl_prof_timer), save :: bpt

    call build(bpt,"compute_nlevs_mg")

    dm = get_dim(ba)

    nlevs_mg = -1

    do i=1,nboxes(ba)

       lo = lwb(get_box(ba,i))
       hi = upb(get_box(ba,i))
       length = hi-lo+1

       do d=1,dm
          temp = length(d)
          rdir = 1
          do while (mod(temp,2) .eq. 0 .and. temp/stag_mg_minwidth .ne. 1)
             temp = temp/2
             rdir = rdir+1
          end do

          if (nlevs_mg .eq. -1) then
             nlevs_mg = rdir
          else          
             nlevs_mg = min(rdir,nlevs_mg)
          end if

       end do

    end do

    call destroy(bpt)

  end subroutine compute_nlevs_mg

end module stag_mg_layout_module
