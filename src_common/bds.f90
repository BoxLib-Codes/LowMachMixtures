module bds_module

  use bl_types
  use multifab_module
  use ml_layout_module
  use define_bc_module
  use convert_stag_module
  use multifab_physbc_extrap_module
  use bc_module
  use probin_common_module, only: advection_type, rhobar

  implicit none

  private

  public :: bds, bds_quad

contains

  subroutine bds(mla,umac,s,s_update,force,s_fc,s_nd,dx,dt,scomp,ncomp, &
                 bc_comp,the_bc_tower,proj_type_in)

    ! s_fc and s_nd are used to set boundary conditions
    ! the input s needs to have ghost cells filled with multifab_physbc_extrap
    ! instead of multifab_physbc

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: umac(:,:)
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: s_update(:)
    type(multifab) , intent(in   ) :: force(:)
    type(multifab) , intent(in   ) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: s_nd(:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    integer        , intent(in   ) :: scomp,ncomp,bc_comp
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ), optional :: proj_type_in

    ! this will hold slx, sly, and slxy
    type(multifab) :: slope(mla%nlevel)

    real(kind=dp_t), pointer ::    sop(:,:,:,:)
    real(kind=dp_t), pointer ::    sup(:,:,:,:)
    real(kind=dp_t), pointer :: slopep(:,:,:,:)
    real(kind=dp_t), pointer ::     fp(:,:,:,:)
    real(kind=dp_t), pointer ::  uadvp(:,:,:,:)
    real(kind=dp_t), pointer ::  vadvp(:,:,:,:)
    real(kind=dp_t), pointer ::  wadvp(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: spz(:,:,:,:)
    real(kind=dp_t), pointer :: sxp(:,:,:,:)
    real(kind=dp_t), pointer :: syp(:,:,:,:)
    real(kind=dp_t), pointer :: szp(:,:,:,:)
    real(kind=dp_t), pointer :: snp(:,:,:,:)

    integer :: dm,ng_s,ng_l,ng_u,ng_v,ng_o,ng_f,ng_e,ng_n,n,i,comp,nlevs,bccomp
    integer :: lo(mla%dim),hi(mla%dim),proj_type

    type(multifab) :: sedge(mla%nlevel,mla%dim)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox, znodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    integer :: zlo(mla%dim), zhi(mla%dim)
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)
    type(box) :: growntilebox
    integer :: gtlo(mla%dim), gthi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"bds")

    ! L2 projection onto EOS?
    ! 0 = do nothing
    ! 1 = assume s holds rho and rho*c
    ! 2 = assume s holds all rho_i's
    if (present(proj_type_in)) then
       proj_type = proj_type_in
    else
       proj_type = 0
    end if

    nlevs = mla%nlevel
    dm = mla%dim

    do n=1,nlevs
       do i=1,dm
          call multifab_build_edge(sedge(n,i),mla%la(n),ncomp,0,i)
       end do
    end do


    if (dm .eq. 2) then
       ! 3 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slxy
       do n=1,nlevs
          call multifab_build(slope(n),mla%la(n),3,1)
       end do
    else if (dm .eq. 3) then
       ! 7 components and 1 ghost cell
       ! component 1 = slx
       ! component 2 = sly
       ! component 3 = slz
       ! component 4 = slxy
       ! component 5 = slxz
       ! component 6 = slyz
       ! component 7 = slxyz
       do n=1,nlevs
          call multifab_build(slope(n),mla%la(n),7,1)
       end do
    end if

    ng_s = s(1)%ng
    ng_u = s_update(1)%ng
    ng_l = slope(1)%ng
    ng_v = umac(1,1)%ng
    ng_o = force(1)%ng
    ng_f = s_fc(1,1)%ng
    ng_e = sedge(1,1)%ng
    ng_n = s_nd(1)%ng

    do comp=scomp,scomp+ncomp-1
       bccomp = bc_comp+comp-scomp

       !$omp parallel private(n,i,mfi,tilebox,tlo,thi,growntilebox,gtlo,gthi) &
       !$omp private(sop,snp,slopep,lo,hi)

       do n=1,nlevs
          call mfiter_build(mfi, s(n), tiling=.true.)

          do while (more_tile(mfi))
             i = get_fab_index(mfi)

             tilebox = get_tilebox(mfi)
             tlo = lwb(tilebox)
             thi = upb(tilebox)
          
             growntilebox = get_growntilebox(mfi,slope(n)%ng)
             gtlo = lwb(growntilebox)
             gthi = upb(growntilebox)

             sop    => dataptr(s(n) , i)
             snp    => dataptr(s_nd(n) , i)
             slopep => dataptr(slope(n), i)
             lo =  lwb(get_box(s_update(n), i))
             hi =  upb(get_box(s_update(n), i))
             select case (dm)
             case (2)
                ! only advancing the tracer
                call bdsslope_2d(lo, hi, &
                                 sop(:,:,1,comp), ng_s, &
                                 slopep(:,:,1,:), ng_l, &
                                 snp(:,:,1,comp), ng_n, &
                                 dx(n,:), &
                                 the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,bccomp), &
                                 tlo,thi,gtlo,gthi) 
             case (3)
                ! only advancing the tracer
                call bdsslope_3d(lo, hi, &
                                 sop(:,:,:,comp), ng_s, &
                                 slopep(:,:,:,:), ng_l, &
                                 snp(:,:,:,comp), ng_n, &
                                 dx(n,:), &
                                 the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,bccomp), &
                                 tlo,thi,gtlo,gthi)

             end select
          end do ! end loop over fabs
       end do ! end loop over levels
       !$omp end parallel


       !$omp parallel private(n,i,mfi,sop,sup,fp,spx,spy,spz,sxp,syp,szp,slopep) &
       !$omp private(uadvp,vadvp,wadvp,lo,hi,tilebox,tlo,thi) &
       !$omp private(xnodalbox,xlo,xhi,ynodalbox,ylo,yhi,znodalbox,zlo,zhi)

       do n=1,nlevs
          call mfiter_build(mfi, s(n), tiling=.true.)


          do while (more_tile(mfi))
             i = get_fab_index(mfi)
             
             tilebox = get_tilebox(mfi)
             tlo = lwb(tilebox)
             thi = upb(tilebox)
             xnodalbox = get_nodaltilebox(mfi,1)
             xlo = lwb(xnodalbox)
             xhi = upb(xnodalbox)
             ynodalbox = get_nodaltilebox(mfi,2)
             ylo = lwb(ynodalbox)
             yhi = upb(ynodalbox)
             znodalbox = get_nodaltilebox(mfi,3)
             zlo = lwb(znodalbox)
             zhi = upb(znodalbox)

             sop    => dataptr(s(n) , i)
             sup    => dataptr(s_update(n), i)
             fp     => dataptr(force(n), i)
             spx => dataptr(s_fc(n,1), i)
             spy => dataptr(s_fc(n,2), i)
             sxp => dataptr(sedge(n,1), i)
             syp => dataptr(sedge(n,2), i)
             slopep => dataptr(slope(n), i)
             uadvp  => dataptr(umac(n,1), i)
             vadvp  => dataptr(umac(n,2), i)
             lo =  lwb(get_box(s_update(n), i))
             hi =  upb(get_box(s_update(n), i))


             select case (dm)
             case (2)
                call bdsconc_2d(lo, hi, &
                                sop(:,:,1,comp), ng_s, sup(:,:,1,comp), ng_u, &
                                fp(:,:,1,comp), ng_o, &
                                slopep(:,:,1,:), ng_l, &
                                uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, &
                                dx(n,:), dt, &
                                spx(:,:,1,comp), spy(:,:,1,comp), ng_f, &
                                sxp(:,:,1,comp-scomp+1), &
                                syp(:,:,1,comp-scomp+1), ng_e, &
                                the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,bccomp), &
                                xlo,xhi,ylo,yhi)
             case (3)
                wadvp  => dataptr(umac(n,3), i)
                spz => dataptr(s_fc(n,3), i)
                szp => dataptr(sedge(n,3), i)
                call bdsconc_3d(lo, hi, &
                                sop(:,:,:,comp), ng_s, sup(:,:,:,comp), ng_u, &
                                fp(:,:,:,comp), ng_o, &
                                slopep(:,:,:,:), ng_l, &
                                uadvp(:,:,:,1), vadvp(:,:,:,1), wadvp(:,:,:,1), ng_v, &
                                dx(n,:), dt, &
                                spx(:,:,:,comp), spy(:,:,:,comp), spz(:,:,:,comp), ng_f, &
                                sxp(:,:,:,comp-scomp+1), &
                                syp(:,:,:,comp-scomp+1), &
                                szp(:,:,:,comp-scomp+1), ng_e, &
                                the_bc_tower%bc_tower_array(n)%adv_bc_level_array(i,:,:,bccomp), &
                                tlo,thi,xlo,xhi,ylo,yhi,zlo,zhi)
             end select
          end do ! end loop over fabs
       end do ! end loop over levels
       !$omp end parallel

    end do ! end loop over components


    ! do L2 projection and increment s_update
    !$omp parallel private(mfi,n,i,xnodalbox,xlo,xhi) &
    !$omp private(ynodalbox,ylo,yhi,znodalbox,zlo,zhi,sxp,syp,szp) &
    !$omp private(uadvp,vadvp,wadvp,lo,hi)
    do n=1,nlevs
       call mfiter_build(mfi, s_update(n), tiling=.true.)
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

          sxp => dataptr(sedge(n,1), i)
          syp => dataptr(sedge(n,2), i)
          uadvp  => dataptr(umac(n,1), i)
          vadvp  => dataptr(umac(n,2), i)
          lo =  lwb(get_box(s_update(n), i))
          hi =  upb(get_box(s_update(n), i))
          select case (dm)
          case (2)
             call bdsupdate_2d(lo, hi, &
                               uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, &
                               sxp(:,:,1,:), syp(:,:,1,:), ng_e, &
                               dx(n,:), ncomp, proj_type,xlo,xhi,ylo,yhi)
          case (3)
             wadvp  => dataptr(umac(n,3), i)
             szp => dataptr(sedge(n,3), i)
             call bdsupdate_3d(lo, hi, &
                               uadvp(:,:,:,1), vadvp(:,:,:,1), wadvp(:,:,:,1), ng_v, &
                               sxp(:,:,:,:), syp(:,:,:,:), szp(:,:,:,:), ng_e, &
                               dx(n,:), ncomp, proj_type,xlo,xhi,ylo,yhi,zlo,zhi)
          end select
       end do ! end loop over fabs
    end do ! end loop over levels
    !$omp end parallel

    !$omp parallel private(mfi,n,i,tilebox,tlo,thi) &
    !$omp private(sup,sxp,syp,szp) &
    !$omp private(uadvp,vadvp,wadvp,lo,hi)
    do n=1,nlevs
       call mfiter_build(mfi, s_update(n), tiling=.true.)
       do while (more_tile(mfi))
          i = get_fab_index(mfi)

          tilebox = get_tilebox(mfi)
          tlo = lwb(tilebox)
          thi = upb(tilebox)

          sup => dataptr(s_update(n), i)
          sxp => dataptr(sedge(n,1), i)
          syp => dataptr(sedge(n,2), i)
          uadvp  => dataptr(umac(n,1), i)
          vadvp  => dataptr(umac(n,2), i)
          lo =  lwb(get_box(s_update(n), i))
          hi =  upb(get_box(s_update(n), i))
          select case (dm)
          case (2)
             call bdsfluxdiv_2d(lo, hi, &
                               sup(:,:,1,scomp:), ng_u, &
                               uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, &
                               sxp(:,:,1,:), syp(:,:,1,:), ng_e, &
                               dx(n,:), ncomp,tlo,thi)
          case (3)
             wadvp  => dataptr(umac(n,3), i)
             szp => dataptr(sedge(n,3), i)
             call bdsfluxdiv_3d(lo, hi, &
                               sup(:,:,:,scomp:), ng_u, &
                               uadvp(:,:,:,1), vadvp(:,:,:,1), wadvp(:,:,:,1), ng_v, &
                               sxp(:,:,:,:), syp(:,:,:,:), szp(:,:,:,:), ng_e, &
                               dx(n,:), ncomp, tlo,thi)
          end select
       end do ! end loop over fabs
    end do ! end loop over levels
    !$omp end parallel

    do n=1,nlevs
       call multifab_destroy(slope(n))
       do i=1,dm
          call multifab_destroy(sedge(n,i))
       end do
    end do

    call destroy(bpt)

  end subroutine bds

  subroutine bdsslope_2d(glo,ghi,s,ng_s,slope,ng_l,s_nd,ng_n,dx,bc,tlo,thi,gtlo,gthi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_s,ng_l,ng_n
    integer        , intent(in   ) :: tlo(:),thi(:),gtlo(:),gthi(:)
    real(kind=dp_t), intent(in   ) ::     s(glo(1)-ng_s:,glo(2)-ng_s:)
    real(kind=dp_t), intent(inout) :: slope(glo(1)-ng_l:,glo(2)-ng_l:,:)
    real(kind=dp_t), intent(in   ) ::  s_nd(glo(1)-ng_n:,glo(2)-ng_n:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:)

    real(kind=dp_t) :: diff(4)
    real(kind=dp_t) :: smin(4)
    real(kind=dp_t) :: smax(4)
    real(kind=dp_t) :: sc(4)

    real(kind=dp_t) :: hx,hy,eps,hxinv,hyinv
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,ll,mm

    ! nodal with one ghost cell
    allocate(sint(tlo(1)-1:thi(1)+2,tlo(2)-1:thi(2)+2))

    hx = dx(1)
    hy = dx(2)

    hxinv = 1/hx
    hyinv = 1/hy

    eps = 1.d-10

    ! bicubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do j = tlo(2)-1,thi(2)+2
       do i = tlo(1)-1,thi(1)+2
          sint(i,j) = (s(i-2,j-2) + s(i-2,j+1) + s(i+1,j-2) + s(i+1,j+1) &
               - 7.d0*(s(i-2,j-1) + s(i-2,j  ) + s(i-1,j-2) + s(i  ,j-2) + & 
                       s(i-1,j+1) + s(i  ,j+1) + s(i+1,j-1) + s(i+1,j  )) &
              + 49.d0*(s(i-1,j-1) + s(i  ,j-1) + s(i-1,j  ) + s(i  ,j  )) ) / 144.d0
       enddo
    enddo

    ! for reservoirs, overwrite boundary nodes and ghost nodes to rho at boundary nodes
    if (tlo(1) .eq. glo(1)) then
    if (bc(1,1) .eq. EXT_DIR) then
       do j=tlo(2)-1,thi(2)+2
          sint(tlo(1)-1:tlo(1),j) = s_nd(tlo(1),j)
       end do
    end if
    end if

    if (thi(1) .eq. ghi(1)) then
    if (bc(1,2) .eq. EXT_DIR) then
       do j=tlo(2)-1,thi(2)+2
          sint(thi(1)+1:thi(1)+2,j) = s_nd(thi(1)+1,j)
       end do
    end if
    end if

    if (tlo(2) .eq. glo(2)) then
    if (bc(2,1) .eq. EXT_DIR) then
       do i=tlo(1)-1,thi(1)+2
          sint(i,tlo(2)-1:tlo(2)) = s_nd(i,tlo(2))
       end do
    end if
    end if

    if (thi(2) .eq. ghi(2)) then
    if (bc(2,2) .eq. EXT_DIR) then
       do i=tlo(1)-1,thi(1)+2
          sint(i,thi(2)+1:thi(2)+2) = s_nd(i,thi(2)+1)
       end do
    end if
    end if

    do j = gtlo(2),gthi(2)
       do i = gtlo(1),gthi(1)

          ! compute initial estimates of slopes from unlimited corner points

          ! sx
          slope(i,j,1) = 0.5d0*(sint(i+1,j+1) + sint(i+1,j  ) - &
                                sint(i  ,j+1) - sint(i  ,j  ) ) * hxinv

          ! sy
          slope(i,j,2) = 0.5d0*(sint(i+1,j+1) - sint(i+1,j  ) + &
                                sint(i  ,j+1) - sint(i  ,j  ) ) * hyinv

          ! sxy
          slope(i,j,3) = ( sint(i+1,j+1) - sint(i+1,j  ) &
                          -sint(i  ,j+1) + sint(i  ,j  ) ) * (hxinv*hyinv)
          
          ! limit slopes
          if (advection_type .eq. 2) then

             ! ++ / sint(i+1,j+1)
             sc(4) = s(i,j) + 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2))  &
                  + 0.25d0*hx*hy*slope(i,j,3)

             ! +- / sint(i+1,j  )
             sc(3) = s(i,j) + 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2))  &
                  - 0.25d0*hx*hy*slope(i,j,3)

             ! -+ / sint(i  ,j+1)
             sc(2) = s(i,j) - 0.5d0*(hx*slope(i,j,1) - hy*slope(i,j,2)) &
                  - 0.25d0*hx*hy*slope(i,j,3)

             ! -- / sint(i  ,j  )
             sc(1) = s(i,j) - 0.5d0*(hx*slope(i,j,1) + hy*slope(i,j,2)) &
                  + 0.25d0*hx*hy*slope(i,j,3)

             ! enforce max/min bounds
             smin(4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smax(4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))

             smin(3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smax(3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))

             smin(2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smax(2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))

             smin(1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             smax(1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

             do mm=1,4
                sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
             enddo

             ! iterative loop
             do ll = 1,3
                sumloc = 0.25d0*(sc(4) + sc(3) + sc(2) + sc(1))
                sumdif = (sumloc - s(i,j))*4.d0
                sgndif = sign(1.d0,sumdif)

                do mm=1,4
                   diff(mm) = (sc(mm) - s(i,j))*sgndif
                enddo

                kdp = 0

                do mm=1,4
                   if (diff(mm) .gt. eps) then
                      kdp = kdp+1
                   end if
                end do

                do mm = 1,4 
                   if (kdp.lt.1) then 
                      div = 1.d0
                   else
                      div = dble(kdp)
                   end if

                   if (diff(mm).gt.eps) then
                      redfac = sumdif*sgndif/div
                      kdp = kdp-1
                   else
                      redfac = 0.d0
                   end if

                   if (sgndif .gt. 0.d0) then
                      redmax = sc(mm) - smin(mm)
                   else
                      redmax = smax(mm) - sc(mm)
                   end if

                   redfac = min(redfac,redmax)
                   sumdif = sumdif - redfac*sgndif
                   sc(mm) = sc(mm) - redfac*sgndif
                enddo
             enddo

             ! final slopes

             ! sx
             slope(i,j,1) = 0.5d0*( sc(4) + sc(3) &
                                   -sc(1) - sc(2)) * hxinv

             ! sy
             slope(i,j,2) = 0.5d0*( sc(4) + sc(2) &
                                   -sc(1) - sc(3)) * hyinv

             ! sxy
             slope(i,j,3) = ( sc(1) + sc(4) &
                             -sc(2) - sc(3) ) * (hxinv * hyinv)

          end if

       enddo
    enddo

    deallocate(sint)

  end subroutine bdsslope_2d

  subroutine bdsslope_3d(glo,ghi,s,ng_s,slope,ng_l,s_nd,ng_n,dx,bc,tlo,thi,gtlo,gthi)

    integer        , intent(in   ) :: glo(:),ghi(:),ng_s,ng_l,ng_n
    integer        , intent(in   ) :: tlo(:),thi(:),gtlo(:),gthi(:)
    real(kind=dp_t), intent(in   ) ::     s(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:)
    real(kind=dp_t), intent(inout) :: slope(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:,:)
    real(kind=dp_t), intent(in   ) ::  s_nd(glo(1)-ng_n:,glo(2)-ng_n:,glo(3)-ng_n:)
    real(kind=dp_t), intent(in   ) :: dx(:)
    integer        , intent(in   ) :: bc(:,:)

    ! local variables
    real(kind=dp_t), allocatable :: sint(:,:,:)

    real(kind=dp_t) :: diff(8)
    real(kind=dp_t) :: smin(8)
    real(kind=dp_t) :: smax(8)
    real(kind=dp_t) :: sc(8)

    real(kind=dp_t) :: c1,c2,c3,c4
    real(kind=dp_t) :: hx,hy,hz,eps,hxinv,hyinv,hzinv
    real(kind=dp_t) :: sumloc,redfac,redmax,div,kdp,sumdif,sgndif
    integer         :: i,j,k,ll,mm

    ! nodal with one ghost cell
    allocate(sint(tlo(1)-1:thi(1)+2,tlo(2)-1:thi(2)+2,tlo(3)-1:thi(3)+2))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    hxinv = 1/hx
    hyinv = 1/hy
    hzinv = 1/hz

    eps = 1.d-10

    c1 = (343.d0/1728.d0)
    c2 = (49.d0 /1728.d0)
    c3 = (7.d0  /1728.d0)
    c4 = (1.d0  /1728.d0)

    ! tricubic interpolation to corner points
    ! (i,j,k) refers to lower corner of cell
    do k = tlo(3)-1,thi(3)+2
       do j = tlo(2)-1,thi(2)+2
          do i = tlo(1)-1,thi(1)+2
             sint(i,j,k) = c1*( s(i  ,j  ,k  ) + s(i-1,j  ,k  ) + s(i  ,j-1,k  ) &
                               +s(i  ,j  ,k-1) + s(i-1,j-1,k  ) + s(i-1,j  ,k-1) &
                               +s(i  ,j-1,k-1) + s(i-1,j-1,k-1) ) &
                          -c2*( s(i-1,j  ,k+1) + s(i  ,j  ,k+1) + s(i-1,j-1,k+1) &
                               +s(i  ,j-1,k+1) + s(i-1,j+1,k  ) + s(i  ,j+1,k  ) &
                               +s(i-2,j  ,k  ) + s(i+1,j  ,k  ) + s(i-2,j-1,k  ) &
                               +s(i+1,j-1,k  ) + s(i-1,j-2,k  ) + s(i  ,j-2,k  ) &
                               +s(i-1,j+1,k-1) + s(i  ,j+1,k-1) + s(i-2,j  ,k-1) &
                               +s(i+1,j  ,k-1) + s(i-2,j-1,k-1) + s(i+1,j-1,k-1) &
                               +s(i-1,j-2,k-1) + s(i  ,j-2,k-1) + s(i-1,j  ,k-2) &
                               +s(i  ,j  ,k-2) + s(i-1,j-1,k-2) + s(i  ,j-1,k-2) ) &
                          +c3*( s(i-1,j+1,k+1) + s(i  ,j+1,k+1) + s(i-2,j  ,k+1) &
                               +s(i+1,j  ,k+1) + s(i-2,j-1,k+1) + s(i+1,j-1,k+1) &
                               +s(i-1,j-2,k+1) + s(i  ,j-2,k+1) + s(i-2,j+1,k  ) &
                               +s(i+1,j+1,k  ) + s(i-2,j-2,k  ) + s(i+1,j-2,k  ) &
                               +s(i-2,j+1,k-1) + s(i+1,j+1,k-1) + s(i-2,j-2,k-1) &
                               +s(i+1,j-2,k-1) + s(i-1,j+1,k-2) + s(i  ,j+1,k-2) &
                               +s(i-2,j  ,k-2) + s(i+1,j  ,k-2) + s(i-2,j-1,k-2) &
                               +s(i+1,j-1,k-2) + s(i-1,j-2,k-2) + s(i  ,j-2,k-2) ) &
                          -c4*( s(i-2,j+1,k+1) + s(i+1,j+1,k+1) + s(i-2,j-2,k+1) &
                               +s(i+1,j-2,k+1) + s(i-2,j+1,k-2) + s(i+1,j+1,k-2) &
                               +s(i-2,j-2,k-2) + s(i+1,j-2,k-2) )
          enddo
       enddo
    enddo

    ! for reservoirs, overwrite boundary nodes and ghost nodes to rho at boundary nodes
    if (tlo(1) .eq. glo(1)) then
    if (bc(1,1) .eq. EXT_DIR) then
       do k=tlo(3)-1,thi(3)+2
       do j=tlo(2)-1,thi(2)+2
          sint(tlo(1)-1:tlo(1),j,k) = s_nd(tlo(1),j,k)
       end do
       end do
    end if
    end if

    if (thi(1) .eq. ghi(1)) then    
    if (bc(1,2) .eq. EXT_DIR) then
       do k=tlo(3)-1,thi(3)+2
       do j=tlo(2)-1,thi(2)+2
          sint(thi(1)+1:thi(1)+2,j,k) = s_nd(thi(1)+1,j,k)
       end do
       end do
    end if
    end if

    if (tlo(2) .eq. glo(2)) then
    if (bc(2,1) .eq. EXT_DIR) then
       do k=tlo(3)-1,thi(3)+2
       do i=tlo(1)-1,thi(1)+2
          sint(i,tlo(2)-1:tlo(2),k) = s_nd(i,tlo(2),k)
       end do
       end do
    end if
    end if

    if (thi(2) .eq. ghi(2)) then
    if (bc(2,2) .eq. EXT_DIR) then
       do k=tlo(3)-1,thi(3)+2
       do i=tlo(1)-1,thi(1)+2
          sint(i,thi(2)+1:thi(2)+2,k) = s_nd(i,thi(2)+1,k)
       end do
       end do
    end if
    end if

    if (tlo(3) .eq. glo(3)) then
    if (bc(3,1) .eq. EXT_DIR) then
       do j=tlo(2)-1,thi(2)+2
       do i=tlo(1)-1,thi(1)+2
          sint(i,j,tlo(3)-1:tlo(3)) = s_nd(i,j,tlo(3))
       end do
       end do
    end if
    end if

    if (thi(3) .eq. ghi(3)) then
    if (bc(3,2) .eq. EXT_DIR) then
       do j=tlo(2)-1,thi(2)+2
       do i=tlo(1)-1,thi(1)+2
          sint(i,j,thi(3)+1:thi(3)+2) = s_nd(i,j,thi(3)+1)
       end do
       end do
    end if
    end if

    do k = gtlo(3),gthi(3)
       do j = gtlo(2),gthi(2)
          do i = gtlo(1),gthi(1) 

             ! compute initial estimates of slopes from unlimited corner points

             ! sx
             slope(i,j,k,1) = 0.25d0*( ( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) * hxinv

             ! sy
             slope(i,j,k,2) = 0.25d0*( ( sint(i  ,j+1,k  ) + sint(i+1,j+1,k  ) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1)) ) * hyinv

             ! sz
             slope(i,j,k,3) = 0.25d0*( ( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                        +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                      -( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                        +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) * hzinv

             ! sxy
             slope(i,j,k,4) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j  ,k+1) &
                                       +sint(i+1,j+1,k  ) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i  ,j+1,k+1)) ) * (hxinv*hyinv)

             ! sxz
             slope(i,j,k,5) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i  ,j+1,k  ) &
                                       +sint(i+1,j  ,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i+1,j  ,k  ) + sint(i+1,j+1,k  ) &
                                       +sint(i  ,j  ,k+1) + sint(i  ,j+1,k+1)) ) * (hxinv*hzinv)

             ! syz
             slope(i,j,k,6) = 0.5d0*( ( sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) &
                                       +sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1)) &
                                     -( sint(i  ,j  ,k+1) + sint(i+1,j  ,k+1) &
                                       +sint(i  ,j+1,k  ) + sint(i+1,j+1,k  )) ) * (hyinv*hzinv)

             ! sxyz
             slope(i,j,k,7) = (-sint(i  ,j  ,k  ) + sint(i+1,j  ,k  ) + sint(i  ,j+1,k  ) &
                               +sint(i  ,j  ,k+1) - sint(i+1,j+1,k  ) - sint(i+1,j  ,k+1) &
                               -sint(i  ,j+1,k+1) + sint(i+1,j+1,k+1) ) * (hxinv*hyinv*hzinv)

             ! limit slopes
             if (advection_type .eq. 2) then

                ! +++ / sint(i+1,j+1,k+1)
                sc(8) = s(i,j,k) &
                     +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                     +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                     +0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! ++- / sint(i+1,j+1,k  )
                sc(7) = s(i,j,k) &
                     +0.5d0*( hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                     +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                     -0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! +-+ / sint(i+1,j  ,k+1)
                sc(6) = s(i,j,k) &
                     +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                     +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                     -0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! +-- / sint(i+1,j  ,k  )
                sc(5) = s(i,j,k) &
                     +0.5d0*( hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                     +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                     +0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! -++ / sint(i  ,j+1,k+1)
                sc(4) = s(i,j,k) &
                     +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                     +0.25d0*(-hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                     -0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! -+- / sint(i  ,j+1,k  )
                sc(3) = s(i,j,k) &
                     +0.5d0*(-hx*slope(i,j,k,1)+hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                     +0.25d0*(-hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                     +0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! --+ / sint(i  ,j  ,k+1)
                sc(2) = s(i,j,k) &
                     +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)+hz*slope(i,j,k,3)) &
                     +0.25d0*( hx*hy*slope(i,j,k,4)-hx*hz*slope(i,j,k,5)-hy*hz*slope(i,j,k,6)) &
                     +0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! ---/ sint(i  ,j  ,k  )
                sc(1) = s(i,j,k) &
                     +0.5d0*(-hx*slope(i,j,k,1)-hy*slope(i,j,k,2)-hz*slope(i,j,k,3)) &
                     +0.25d0*( hx*hy*slope(i,j,k,4)+hx*hz*slope(i,j,k,5)+hy*hz*slope(i,j,k,6)) &
                     -0.125d0*hx*hy*hz*slope(i,j,k,7)

                ! enforce max/min bounds
                smin(8) = min(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                     s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))
                smax(8) = max(s(i,j,k),s(i+1,j,k),s(i,j+1,k),s(i,j,k+1), &
                     s(i+1,j+1,k),s(i+1,j,k+1),s(i,j+1,k+1),s(i+1,j+1,k+1))

                smin(7) = min(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                     s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))
                smax(7) = max(s(i,j,k-1),s(i+1,j,k-1),s(i,j+1,k-1),s(i,j,k), &
                     s(i+1,j+1,k-1),s(i+1,j,k),s(i,j+1,k),s(i+1,j+1,k))

                smin(6) = min(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                     s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))
                smax(6) = max(s(i,j-1,k),s(i+1,j-1,k),s(i,j,k),s(i,j-1,k+1), &
                     s(i+1,j,k),s(i+1,j-1,k+1),s(i,j,k+1),s(i+1,j,k+1))

                smin(5) = min(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                     s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))
                smax(5) = max(s(i,j-1,k-1),s(i+1,j-1,k-1),s(i,j,k-1),s(i,j-1,k), &
                     s(i+1,j,k-1),s(i+1,j-1,k),s(i,j,k),s(i+1,j,k))

                smin(4) = min(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                     s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))
                smax(4) = max(s(i-1,j,k),s(i,j,k),s(i-1,j+1,k),s(i-1,j,k+1), &
                     s(i,j+1,k),s(i,j,k+1),s(i-1,j+1,k+1),s(i,j+1,k+1))

                smin(3) = min(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                     s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))
                smax(3) = max(s(i-1,j,k-1),s(i,j,k-1),s(i-1,j+1,k-1),s(i-1,j,k), &
                     s(i,j+1,k-1),s(i,j,k),s(i-1,j+1,k),s(i,j+1,k))

                smin(2) = min(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                     s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))
                smax(2) = max(s(i-1,j-1,k),s(i,j-1,k),s(i-1,j,k),s(i-1,j-1,k+1), &
                     s(i,j,k),s(i,j-1,k+1),s(i-1,j,k+1),s(i,j,k+1))

                smin(1) = min(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                     s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))
                smax(1) = max(s(i-1,j-1,k-1),s(i,j-1,k-1),s(i-1,j,k-1),s(i-1,j-1,k), &
                     s(i,j,k-1),s(i,j-1,k),s(i-1,j,k),s(i,j,k))

                do mm=1,8
                   sc(mm) = max(min(sc(mm), smax(mm)), smin(mm))
                enddo

                ! iterative loop
                do ll = 1,4
                   sumloc = 0.125d0*(sc(1)+sc(2)+sc(3)+sc(4)+sc(5)+sc(6)+sc(7)+sc(8))
                   sumdif = (sumloc - s(i,j,k))*8.d0
                   sgndif = sign(1.d0,sumdif)

                   do mm=1,8
                      diff(mm) = (sc(mm) - s(i,j,k))*sgndif
                   enddo

                   kdp = 0

                   do mm=1,8
                      if (diff(mm) .gt. eps) then
                         kdp = kdp+1
                      end if
                   end do

                   do mm = 1,8
                      if (kdp.lt.1) then 
                         div = 1.d0
                      else
                         div = dble(kdp)
                      end if

                      if (diff(mm).gt.eps) then
                         redfac = sumdif*sgndif/div
                         kdp = kdp-1
                      else
                         redfac = 0.d0
                      end if

                      if (sgndif .gt. 0.d0) then
                         redmax = sc(mm) - smin(mm)
                      else
                         redmax = smax(mm) - sc(mm)
                      end if

                      redfac = min(redfac,redmax)
                      sumdif = sumdif - redfac*sgndif
                      sc(mm) = sc(mm) - redfac*sgndif
                   enddo
                enddo

                ! final slopes

                ! sx
                slope(i,j,k,1) = 0.25d0*( ( sc(5) + sc(7) &
                                           +sc(6) + sc(8)) &
                                         -( sc(1) + sc(3) &
                                           +sc(2) + sc(4)) ) * hxinv

                ! sy
                slope(i,j,k,2) = 0.25d0*( ( sc(3) + sc(7) &
                                           +sc(4) + sc(8)) &
                                         -( sc(1) + sc(5) &
                                           +sc(2) + sc(6)) ) * hyinv

                ! sz
                slope(i,j,k,3) = 0.25d0*( ( sc(2) + sc(6) &
                                           +sc(4) + sc(8)) &
                                         -( sc(1) + sc(5) &
                                           +sc(3) + sc(7)) ) * hzinv

                ! sxy
                slope(i,j,k,4) = 0.5d0*( ( sc(1) + sc(2) &
                                          +sc(7) + sc(8)) &
                                        -( sc(5) + sc(6) &
                                          +sc(3) + sc(4)) ) * (hxinv*hyinv)

                ! sxz
                slope(i,j,k,5) = 0.5d0*( ( sc(1) + sc(3) &
                                          +sc(6) + sc(8)) &
                                        -( sc(5) + sc(7) &
                                          +sc(2) + sc(4)) ) * (hxinv*hzinv)

                ! syz
                slope(i,j,k,6) = 0.5d0*( ( sc(1) + sc(5) &
                                          +sc(4) + sc(8)) &
                                        -( sc(2) + sc(6) &
                                          +sc(3) + sc(7)) ) * (hyinv*hzinv)

                ! sxyz
                slope(i,j,k,7) = (-sc(1) + sc(5) + sc(3) &
                                          +sc(2) - sc(7) - sc(6) &
                                          -sc(4) + sc(8) ) * (hxinv*hyinv*hzinv)

             endif

          enddo
       enddo
    enddo

  end subroutine bdsslope_3d

  subroutine bdsconc_2d(glo,ghi,s,ng_s,s_update,ng_u,force,ng_o, &
                        slope,ng_l,uadv,vadv,ng_v,dx,dt,sx,sy,ng_f, &
                        sedgex,sedgey,ng_e,bc,xlo,xhi,ylo,yhi)

    integer        ,intent(in   ) :: glo(:),ghi(:),ng_s,ng_u,ng_l,ng_v,ng_o,ng_f,ng_e
    integer        ,intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
    real(kind=dp_t),intent(in   ) ::        s(glo(1)-ng_s:,glo(2)-ng_s:)
    real(kind=dp_t),intent(inout) :: s_update(glo(1)-ng_u:,glo(2)-ng_u:)
    real(kind=dp_t),intent(in   ) ::    force(glo(1)-ng_o:,glo(2)-ng_o:)
    real(kind=dp_t),intent(in   ) ::    slope(glo(1)-ng_l:,glo(2)-ng_l:,:)
    real(kind=dp_t),intent(in   ) ::     uadv(glo(1)-ng_v:,glo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(glo(1)-ng_v:,glo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) ::       sx(glo(1)-ng_f:,glo(2)-ng_f:)
    real(kind=dp_t),intent(in   ) ::       sy(glo(1)-ng_f:,glo(2)-ng_f:)
    real(kind=dp_t),intent(inout) ::   sedgex(glo(1)-ng_e:,glo(2)-ng_e:)
    real(kind=dp_t),intent(inout) ::   sedgey(glo(1)-ng_e:,glo(2)-ng_e:)
    real(kind=dp_t),intent(in   ) :: dx(:),dt
    integer        ,intent(in   ) :: bc(:,:)

    ! local variables
    real(kind=dp_t) :: hx,hy,hxs,hys,gamp,gamm,hxinv,hyinv
    real(kind=dp_t) :: vtrans,stem,vaddif,vdif
    real(kind=dp_t) :: isign, jsign
    real(kind=dp_t) :: u1,u2,v1,v2,uu,vv

    integer i,j,iup,jup

    hx = dx(1)
    hy = dx(2)

    hxinv = 1/hx
    hyinv = 1/hy

    do j = xlo(2),xhi(2) 
       do i = xlo(1)-1,xhi(1)-1

          ! ******************************* 
          ! calculate Gamma plus for flux F

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j+1)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j
             jsign = 1.d0
             u2 = uadv(i+1,j)
          else 
             jup   = j+1
             jsign = -1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j+1) .gt. 0) then
                u2 = uadv(i+1,j+1)
             end if
          end if

          vv = vadv(iup,j+1)

          hxs = hx*isign
          hys = hy*jsign

          gamp = s(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slope(iup,jup,1) +   &
               (hys*.5 -    vv*dt/3.d0)*slope(iup,jup,2) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma plus for flux F
          ! ****************************************

          ! *****************************************
          ! calculate Gamma minus for flux F

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j-1
             jsign = 1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j-1) .gt. 0) then
                u2 = uadv(i+1,j-1)
             end if
          else 
             jup   = j
             jsign = -1.d0
             u2 = uadv(i+1,j)
          end if

          vv = vadv(iup,j)

          hxs = hx*isign
          hys = hy*jsign

          gamm = s(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slope(iup,jup,1) +   &
               (hys*.5 -    vv*dt/3.d0)*slope(iup,jup,2) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma minus for flux F
          ! ****************************************

          ! *********************************
          ! calculate sedgex

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vdif = 0.5d0*dt*(vadv(iup,j+1)*gamp -  &
               vadv(iup,j)*gamm ) * hyinv
          stem = s(iup,j) + (isign*hx - uadv(i+1,j)*dt)*0.5d0*slope(iup,j,1)
          vaddif = stem*0.5d0*dt*( &
               uadv(iup+1,j) - uadv(iup,j))*hxinv
          sedgex(i+1,j) = stem - vdif - vaddif + 0.5d0*dt*force(iup,j)

          ! end of calculation of sedgex
          ! *************************************

       enddo
    enddo

    if (xlo(1) .eq. glo(1)) then
    if (bc(1,1) .eq. EXT_DIR) then
       sedgex(xlo(1),xlo(2):xhi(2)) = sx(xlo(1),xlo(2):xhi(2))
    end if
    end if

    if (xhi(1) .eq. ghi(1)+1) then
    if (bc(1,2) .eq. EXT_DIR) then
       sedgex(xhi(1),xlo(2):xhi(2)) = sx(xhi(1),xlo(2):xhi(2))
    end if
    end if

    do j = ylo(2)-1,yhi(2)-1
       do i = ylo(1),yhi(1)

          ! ********************************** 
          ! calculate Gamma plus for flux G

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i+1,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i
             isign = 1.d0
             v2 = vadv(i,j+1)
          else
             iup   = i+1
             isign = -1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i+1,j+1) .gt. 0) then
                v2 = vadv(i+1,j+1)
             end if
          end if

          uu = uadv(i+1,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamp = s(iup,jup)+ &
               (hys*.5 - (v1+v2)*dt/3.)*slope(iup,jup,2) +   &
               (hxs*.5 - uu*dt/3.)*slope(iup,jup,1) + &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma plus for flux G
          ! ****************************************

          ! *****************************************
          ! calculate Gamma minus for flux G

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i-1
             isign = 1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i-1,j+1) .gt. 0) then
                v2 = vadv(i-1,j+1)
             end if
          else
             iup   = i
             isign = -1.d0
             v2 = vadv(i,j+1)
          end if

          uu = uadv(i,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamm = s(iup,jup) +    &
               (hys*.5 - (v1+v2)*dt/3.)*slope(iup,jup,2) +    &
               (hxs*.5 - uu*dt/3.)*slope(iup,jup,1) +   &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slope(iup,jup,3)/12.d0

          ! end of calculation of Gamma minus for flux G
          ! ****************************************

          ! *********************************
          ! calculate sedgey

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vdif = 0.5d0*dt* &
               (uadv(i+1,jup)*gamp-uadv(i,jup)*gamm)*hxinv
          stem = s(i,jup) + (jsign*hy - vadv(i,j+1)*dt)*0.5d0*slope(i,jup,2)
          vaddif = stem*0.5d0*dt*(vadv(i,jup+1) - vadv(i,jup))*hyinv
          sedgey(i,j+1) = stem - vdif - vaddif + 0.5d0*dt*force(i,jup)

          ! end of calculation of sedgey
          ! *************************************

       enddo
    enddo

    if (ylo(2) .eq. glo(2)) then
    if (bc(2,1) .eq. EXT_DIR) then
       sedgey(ylo(1):yhi(1),ylo(2)) = sy(ylo(1):yhi(1),ylo(2))
    end if
    end if

    if (yhi(2) .eq. ghi(2)+1) then
    if (bc(2,2) .eq. EXT_DIR) then
       sedgey(ylo(1):yhi(1),yhi(2)) = sy(ylo(1):yhi(1),yhi(2))
    end if
    end if

  end subroutine bdsconc_2d

  subroutine bdsupdate_2d(lo,hi,uadv,vadv,ng_v, &
                          sedgex,sedgey,ng_e,dx,ncomp,proj_type,xlo,xhi,ylo,yhi)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_v,ng_e,ncomp,proj_type
    integer        ,intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:)
    real(kind=dp_t),intent(in   ) ::     uadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(inout) ::   sedgex(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgey(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t),intent(in   ) :: dx(:)

    ! local variables
    integer comp,i,j

    real(kind=dp_t) :: w(ncomp), rhobar_sq, delta_eos, temp


    if (proj_type .eq. 1 .or. proj_type .eq. 2) then

       ! L2 projection: x-faces
       do j = xlo(2),xhi(2) 
       do i = xlo(1),xhi(1) 
          
          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains rho_i instead of (rho,rho_1)
             temp = sedgex(i,j,1)
             sedgex(i,j,1) = sedgex(i,j,2)
             sedgex(i,j,2) = temp - sedgex(i,j,1)
          end if

          ! compute mass fractions, w_i = rho_i / rho
          temp = 0.d0
          do comp=1,ncomp
             temp = temp + sedgex(i,j,comp)
          end do
          do comp=1,ncomp
             w(comp) = sedgex(i,j,comp)/temp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do comp=1,ncomp
             rhobar_sq = rhobar_sq + w(comp)/rhobar(comp)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do comp=1,ncomp
             delta_eos = delta_eos + sedgex(i,j,comp)/rhobar(comp)
          end do

          do comp=1,ncomp
             sedgex(i,j,comp) = sedgex(i,j,comp) - w(comp)*(rhobar_sq/rhobar(comp))*delta_eos
          end do

          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains (rho,rho_1) instead of rho_i
             temp = sedgex(i,j,1)
             sedgex(i,j,1) = sedgex(i,j,1) + sedgex(i,j,2)
             sedgex(i,j,2) = temp
          end if
          
       enddo
       enddo

       ! L2 projection: y-faces
       do j = ylo(2),yhi(2) 
       do i = ylo(1),yhi(1) 
          
          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains rho_i instead of (rho,rho_1)
             temp = sedgey(i,j,1)
             sedgey(i,j,1) = sedgey(i,j,2)
             sedgey(i,j,2) = temp - sedgey(i,j,1)
          end if

          ! compute mass fractions, w_i = rho_i / rho
          temp = 0.d0
          do comp=1,ncomp
             temp = temp + sedgey(i,j,comp)
          end do
          do comp=1,ncomp
             w(comp) = sedgey(i,j,comp)/temp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do comp=1,ncomp
             rhobar_sq = rhobar_sq + w(comp)/rhobar(comp)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do comp=1,ncomp
             delta_eos = delta_eos + sedgey(i,j,comp)/rhobar(comp)
          end do

          do comp=1,ncomp
             sedgey(i,j,comp) = sedgey(i,j,comp) - w(comp)*(rhobar_sq/rhobar(comp))*delta_eos
          end do

          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains (rho,rho_1) instead of rho_i
             temp = sedgey(i,j,1)
             sedgey(i,j,1) = sedgey(i,j,1) + sedgey(i,j,2)
             sedgey(i,j,2) = temp
          end if

       enddo
       enddo

    end if



  end subroutine bdsupdate_2d

  subroutine bdsfluxdiv_2d(lo,hi,s_update,ng_u,uadv,vadv,ng_v, &
                          sedgex,sedgey,ng_e,dx,ncomp,tlo,thi)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_u,ng_v,ng_e,ncomp,tlo(:),thi(:)
    real(kind=dp_t),intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,:)
    real(kind=dp_t),intent(in   ) ::     uadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(inout) ::   sedgex(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgey(lo(1)-ng_e:,lo(2)-ng_e:,:)
    real(kind=dp_t),intent(in   ) :: dx(:)

    ! local variables
    integer comp,i,j

    real(kind=dp_t) :: w(ncomp), rhobar_sq, delta_eos, temp, dxinv

    dxinv = 1.d0/dx(1) !assumed all dx equal

    ! advance solution
    ! conservative update

    do comp=1,ncomp
       do j = tlo(2),thi(2) 
       do i = tlo(1),thi(1) 
          s_update(i,j,comp) = s_update(i,j,comp) - (  &
               (sedgex(i+1,j,comp)*uadv(i+1,j)-sedgex(i,j,comp)*uadv(i,j))*dxinv +  &
               (sedgey(i,j+1,comp)*vadv(i,j+1)-sedgey(i,j,comp)*vadv(i,j))*dxinv )
       enddo
       enddo
    enddo
  end subroutine bdsfluxdiv_2d


  subroutine bdsconc_3d(glo,ghi,s,ng_s,s_update,ng_u,force,ng_o, &
                        slope,ng_l,uadv,vadv,wadv,ng_v,dx,dt,sx,sy,sz,ng_f, &
                        sedgex,sedgey,sedgez,ng_e,bc,tlo,thi, &
                        xlo,xhi,ylo,yhi,zlo,zhi) 

    integer        ,intent(in   ) :: glo(:),ghi(:),ng_s,ng_l,ng_v,ng_u,ng_o,ng_f,ng_e
    integer        ,intent(in   ) :: tlo(:),thi(:),xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t),intent(in   ) ::        s(glo(1)-ng_s:,glo(2)-ng_s:,glo(3)-ng_s:)
    real(kind=dp_t),intent(inout) :: s_update(glo(1)-ng_u:,glo(2)-ng_u:,glo(3)-ng_u:)
    real(kind=dp_t),intent(in   ) ::    force(glo(1)-ng_o:,glo(2)-ng_o:,glo(3)-ng_o:)
    real(kind=dp_t),intent(in   ) ::    slope(glo(1)-ng_l:,glo(2)-ng_l:,glo(3)-ng_l:,:)
    real(kind=dp_t),intent(in   ) ::     uadv(glo(1)-ng_v:,glo(2)-ng_v:,glo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(glo(1)-ng_v:,glo(2)-ng_v:,glo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     wadv(glo(1)-ng_v:,glo(2)-ng_v:,glo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::       sx(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:)
    real(kind=dp_t),intent(in   ) ::       sy(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:)
    real(kind=dp_t),intent(in   ) ::       sz(glo(1)-ng_f:,glo(2)-ng_f:,glo(3)-ng_f:)
    real(kind=dp_t),intent(inout) ::   sedgex(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t),intent(inout) ::   sedgey(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t),intent(inout) ::   sedgez(glo(1)-ng_e:,glo(2)-ng_e:,glo(3)-ng_e:)
    real(kind=dp_t),intent(in   ) :: dx(:),dt
    integer        ,intent(in   ) :: bc(:,:)

    ! local variables
    integer i,j,k,ioff,joff,koff,ll

    real(kind=dp_t), allocatable :: ux(:,:,:)
    real(kind=dp_t), allocatable :: vy(:,:,:)
    real(kind=dp_t), allocatable :: wz(:,:,:)

    real(kind=dp_t) :: isign,jsign,ksign,hx,hy,hz,hxinv,hyinv,hzinv
    real(kind=dp_t) :: del(3),p1(3),p2(3),p3(3),p4(3)
    real(kind=dp_t) :: val1,val2,val3,val4,val5
    real(kind=dp_t) :: u,v,w,uu,vv,ww,gamma,gamma2
    real(kind=dp_t) :: dt2,dt3,dt4,half,sixth

    allocate(ux(tlo(1)-1:thi(1)+1,tlo(2)-1:thi(2)+1,tlo(3)-1:thi(3)+1))
    allocate(vy(tlo(1)-1:thi(1)+1,tlo(2)-1:thi(2)+1,tlo(3)-1:thi(3)+1))
    allocate(wz(tlo(1)-1:thi(1)+1,tlo(2)-1:thi(2)+1,tlo(3)-1:thi(3)+1))

    hx = dx(1)
    hy = dx(2)
    hz = dx(3)

    hxinv = 1/hx
    hyinv = 1/hy
    hzinv = 1/hz

    dt2 = dt/2.d0
    dt3 = dt/3.d0
    dt4 = dt/4.d0

    half = 0.5d0
    sixth = 1.d0/6.d0

    ! compute cell-centered ux, vy, and wz
    do k=tlo(3)-1,thi(3)+1
       do j=tlo(2)-1,thi(2)+1
          do i=tlo(1)-1,thi(1)+1
             ux(i,j,k) = (uadv(i+1,j,k) - uadv(i,j,k)) * hxinv
             vy(i,j,k) = (vadv(i,j+1,k) - vadv(i,j,k)) * hyinv
             wz(i,j,k) = (wadv(i,j,k+1) - wadv(i,j,k)) * hzinv
          end do
       end do
    end do

    ! compute sedgex on x-faces
    do k=xlo(3),xhi(3)
       do j=xlo(2),xhi(2)
          do i=xlo(1),xhi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgex without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ! centroid of rectangular volume
             del(1) = isign*0.5d0*hx - 0.5d0*uadv(i,j,k)*dt
             del(2) = 0.d0
             del(3) = 0.d0
             call eval(s(i+ioff,j,k),slope(i+ioff,j,k,:),del,sedgex(i,j,k))
    
             ! source term
             sedgex(i,j,k) = sedgex(i,j,k) - dt2*sedgex(i,j,k)*ux(i+ioff,j,k) &
                  + dt2*force(i+ioff,j,k)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k) .gt. 0) then
                u = uadv(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j+1,k)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i+ioff,j+1,k)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k) .gt. 0) then
                u = uadv(i,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j+joff,k) + gamma*vy(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             vv = 0.d0
             if (vadv(i+ioff,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j,k+koff) .gt. 0) then
                u = uadv(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k+1)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i+ioff,j,k+1)
             sedgex(i,j,k) = sedgex(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             u = 0.d0
             if (uadv(i,j,k)*uadv(i,j,k+koff) .gt. 0) then
                u = uadv(i,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz

             p3(1) = isign*0.5d0*hx - u*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*ux(i+ioff,j,k+koff) + gamma*wz(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             uu = 0.d0
             if (uadv(i,j,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             ww = 0.d0
             if (wadv(i+ioff,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i+ioff,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgex with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i+ioff,j,k)
             sedgex(i,j,k) = sedgex(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo

    if (xlo(1) .eq. glo(1)) then
    if (bc(1,1) .eq. EXT_DIR) then
       sedgex(xlo(1),xlo(2):xhi(2),xlo(3):xhi(3)) = sx(xlo(1),xlo(2):xhi(2),xlo(3):xhi(3))
    end if
    end if

    if (xhi(1) .eq. ghi(1)+1) then
    if (bc(1,2) .eq. EXT_DIR) then
       sedgex(xhi(1),xlo(2):xhi(2),xlo(3):xhi(3)) = sx(xhi(1),xlo(2):xhi(2),xlo(3):xhi(3))
    end if
    end if


    ! compute sedgey on y-faces    
    do k=ylo(3),yhi(3)
       do j=ylo(2),yhi(2)
          do i=ylo(1),yhi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgey without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (vadv(i,j,k) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             del(1) = 0.d0
             del(2) = jsign*0.5d0*hy - 0.5d0*vadv(i,j,k)*dt
             del(3) = 0.d0
             call eval(s(i,j+joff,k),slope(i,j+joff,k,:),del,sedgey(i,j,k))

             ! source term
             sedgey(i,j,k) = sedgey(i,j,k) - dt2*sedgey(i,j,k)*vy(i,j+joff,k) &
                  + dt2*force(i,j+joff,k)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k) .gt. 0) then
                v = vadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i+1,j+joff,k)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * uadv(i+1,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k) .gt. 0) then
                v = vadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = 0.d0

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = 0.d0

             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = 0.d0

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j+joff,k),slope(i+ioff,j+joff,k,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i+ioff,j+joff,k) + gamma*ux(i+ioff,j+joff,k))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k+1)

             gamma = gamma - dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i+ioff,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             uu = 0.d0
             if (uadv(i,j+joff,k)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - wadv(i+ioff,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * wadv(i+ioff,j+joff,k)

             gamma = gamma + dt*gamma2/(3.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i,j+joff,k+1) .gt. 0) then
                ksign = 1.d0
                koff = 0
             else
                ksign = -1.d0
                koff = 1
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i,j,k+koff) .gt. 0) then
                v = vadv(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z+} with \Gamma^{z+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k+1)*wadv(i+ioff,j+joff,k+1) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k+1)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k+1)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i,j+joff,k+1)
             sedgey(i,j,k) = sedgey(i,j,k) - dt*gamma/(2.d0*hz)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{z-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (wadv(i,j+joff,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             v = 0.d0
             if (vadv(i,j,k)*vadv(i,j,k+koff) .gt. 0) then
                v = vadv(i,j,k+koff)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - v*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*vy(i,j+joff,k+koff) + gamma*wz(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{z-} with \Gamma^{z-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             vv = 0.d0
             if (vadv(i,j,k)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             ww = 0.d0
             if (wadv(i,j+joff,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p2(3) = ksign*0.5d0*hz
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j+joff,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgey with \Gamma^{z-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * wadv(i,j+joff,k)
             sedgey(i,j,k) = sedgey(i,j,k) + dt*gamma/(2.d0*hz)

          enddo
       enddo
    enddo


    if (ylo(2) .eq. glo(2)) then    
    if (bc(2,1) .eq. EXT_DIR) then
       sedgey(ylo(1):yhi(1),ylo(2),ylo(3):yhi(3)) = sy(ylo(1):yhi(1),ylo(2),ylo(3):yhi(3))
    end if
    end if

    if (yhi(2) .eq. ghi(2)+1) then
    if (bc(2,2) .eq. EXT_DIR) then
       sedgey(ylo(1):yhi(1),yhi(2),ylo(3):yhi(3)) = sy(ylo(1):yhi(1),yhi(2),ylo(3):yhi(3))
    end if
    end if


    ! compute sedgez on z-faces
    do k=zlo(3),zhi(3)
       do j=zlo(2),zhi(2)
          do i=zlo(1),zhi(1)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute sedgez without transverse corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             ! centroid of rectangular volume
             if (wadv(i,j,k) .gt. 0) then
                ksign = 1.d0
                koff = -1
             else
                ksign = -1.d0
                koff = 0
             endif

             del(1) = 0.d0
             del(2) = 0.d0
             del(3) = ksign*0.5d0*hz - 0.5d0*wadv(i,j,k)*dt
             call eval(s(i,j,k+koff),slope(i,j,k+koff,:),del,sedgez(i,j,k))

             ! source term
             sedgez(i,j,k) = sedgez(i,j,k) - dt2*sedgez(i,j,k)*wz(i,j,k+koff) &
                  + dt2*force(i,j,k+koff)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j,k) .gt. 0) then
                w = wadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - uadv(i+1,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x+} with \Gamma^{x+,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i+1,j,k+koff)*uadv(i+1,j+joff,k+koff) .gt. 0) then
                uu = uadv(i+1,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i+1,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{x-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j,k) .gt. 0) then
                w = wadv(i+ioff,j,k)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = 0.d0
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = 0.d0
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = isign*0.5d0*hx - uadv(i,j,k+koff)*dt
             p3(2) = 0.d0
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i+ioff,j,k+koff),slope(i+ioff,j,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i+ioff,j,k+koff) + gamma*ux(i+ioff,j,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j+1,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{x-} with \Gamma^{x-,y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i+ioff,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             uu = 0.d0
             if (uadv(i,j,k+koff)*uadv(i,j+joff,k+koff) .gt. 0) then
                uu = uadv(i,j+joff,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx - uadv(i,j+joff,k)*dt
             p3(2) = jsign*0.5d0*hy
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uu*dt
             p4(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k+koff)*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * vadv(i+ioff,j,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * uadv(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y+} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i,j+1,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = 0
             else
                jsign = -1.d0
                joff = 1
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i,j+joff,k) .gt. 0) then
                w = wadv(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vadv(i,j+1,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y+} with \Gamma^{y+,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j+1,k+koff)*vadv(i+ioff,j+1,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j+1,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j+1,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             
             gamma = gamma * vadv(i,j+1,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) - dt*gamma/(2.d0*hy)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! compute \Gamma^{y-} without corner corrections
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (vadv(i,j,k+koff) .gt. 0) then
                jsign = 1.d0
                joff = -1
             else
                jsign = -1.d0
                joff = 0
             endif

             w = 0.d0
             if (wadv(i,j,k)*wadv(i,j+joff,k) .gt. 0) then
                w = wadv(i,j+joff,k)
             endif

             p1(1) = 0.d0
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = 0.d0
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p3(1) = 0.d0
             p3(2) = jsign*0.5d0*hy - vadv(i,j,k+koff)*dt
             p3(3) = ksign*0.5d0*hz - w*dt

             do ll=1,3
                del(ll) = (p2(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = (p1(ll)+p3(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll))/2.d0
             end do
             call eval(s(i,j+joff,k+koff),slope(i,j+joff,k+koff,:),del,val3)

             ! average these centroid values to get the average value
             gamma = (val1+val2+val3)/3.d0

             ! source term
             gamma = gamma - dt3*(gamma*wz(i,j+joff,k+koff) + gamma*vy(i,j+joff,k+koff))

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x+}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i+1,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = 0
             else
                isign = -1.d0
                ioff = 1
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i+1,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i+1,j+joff,k+koff)

             gamma = gamma - dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct \Gamma^{y-} with \Gamma^{y-,x-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             if (uadv(i,j+joff,k+koff) .gt. 0) then
                isign = 1.d0
                ioff = -1
             else
                isign = -1.d0
                ioff = 0
             endif

             ww = 0.d0
             if (wadv(i,j,k)*wadv(i+ioff,j+joff,k) .gt. 0) then
                ww = wadv(i+ioff,j+joff,k)
             endif

             vv = 0.d0
             if (vadv(i,j,k+koff)*vadv(i+ioff,j,k+koff) .gt. 0) then
                vv = vadv(i+ioff,j,k+koff)
             endif

             p1(1) = isign*0.5d0*hx
             p1(2) = jsign*0.5d0*hy
             p1(3) = ksign*0.5d0*hz

             p2(1) = isign*0.5d0*hx
             p2(2) = jsign*0.5d0*hy
             p2(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt
             
             p3(1) = isign*0.5d0*hx
             p3(2) = jsign*0.5d0*hy - vadv(i+ioff,j,k)*dt
             p3(3) = ksign*0.5d0*hz - wadv(i,j,k)*dt

             p4(1) = isign*0.5d0*hx - uadv(i,j+joff,k+koff)*dt
             p4(2) = jsign*0.5d0*hy - vv*dt
             p4(3) = ksign*0.5d0*hz - ww*dt

             do ll=1,3
                del(ll) = (p1(ll)+p2(ll)+p3(ll)+p4(ll))/4.d0
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val1)

             do ll=1,3
                del(ll) = half*p1(ll) + sixth*(p2(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val2)

             do ll=1,3
                del(ll) = half*p2(ll) + sixth*(p1(ll)+p3(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val3)

             do ll=1,3
                del(ll) = half*p3(ll) + sixth*(p2(ll)+p1(ll)+p4(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val4)

             do ll=1,3
                del(ll) = half*p4(ll) + sixth*(p2(ll)+p3(ll)+p1(ll))
             end do
             call eval(s(i+ioff,j+joff,k+koff),slope(i+ioff,j+joff,k+koff,:),del,val5)

             gamma2 = -0.8d0*val1 + 0.45d0*(val2+val3+val4+val5)

             ! source term
             gamma2 = gamma2 - dt4 * ( gamma2*ux(i+ioff,j+joff,k+koff) &
                                      +gamma2*vy(i+ioff,j+joff,k+koff) &
                                      +gamma2*wz(i+ioff,j+joff,k+koff))

             gamma2 = gamma2 * uadv(i,j+joff,k+koff)

             gamma = gamma + dt*gamma2/(3.d0*hx)

             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             ! correct sedgez with \Gamma^{y-}
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

             gamma = gamma * vadv(i,j,k+koff)
             sedgez(i,j,k) = sedgez(i,j,k) + dt*gamma/(2.d0*hy)

          enddo
       enddo
    enddo

    if (zlo(3) .eq. glo(3)) then
    if (bc(3,1) .eq. EXT_DIR) then
       sedgez(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3)) = sz(zlo(1):zhi(1),zlo(2):zhi(2),zlo(3))
    end if
    end if

    if (zhi(3) .eq. ghi(3)+1) then
    if (bc(3,2) .eq. EXT_DIR) then
       sedgez(zlo(1):zhi(1),zlo(2):zhi(2),zhi(3)) = sz(zlo(1):zhi(1),zlo(2):zhi(2),zhi(3))
    end if
    end if


    deallocate(ux,vy,wz)




  end subroutine bdsconc_3d

  subroutine bdsupdate_3d(lo,hi,uadv,vadv,wadv,ng_v, &
                          sedgex,sedgey,sedgez,ng_e,dx,ncomp,proj_type, &
                          xlo,xhi,ylo,yhi,zlo,zhi)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_v,ng_e,ncomp,proj_type
    integer        ,intent(in   ) :: xlo(:),xhi(:),ylo(:),yhi(:),zlo(:),zhi(:)
    real(kind=dp_t),intent(in   ) ::     uadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     wadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(inout) ::   sedgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(in   ) :: dx(:)

    ! local variables
    integer comp,i,j,k

    real(kind=dp_t) :: w(ncomp), rhobar_sq, delta_eos, temp

    if (proj_type .eq. 1 .or. proj_type .eq. 2) then

       ! L2 projection: x-faces
       do k = xlo(3),xhi(3)
       do j = xlo(2),xhi(2) 
       do i = xlo(1),xhi(1) 
          
          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains rho_i instead of (rho,rho_1)
             temp = sedgex(i,j,k,1)
             sedgex(i,j,k,1) = sedgex(i,j,k,2)
             sedgex(i,j,k,2) = temp - sedgex(i,j,k,1)
          end if

          ! compute mass fractions, w_i = rho_i / rho
          temp = 0.d0
          do comp=1,ncomp
             temp = temp + sedgex(i,j,k,comp)
          end do
          do comp=1,ncomp
             w(comp) = sedgex(i,j,k,comp)/temp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do comp=1,ncomp
             rhobar_sq = rhobar_sq + w(comp)/rhobar(comp)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do comp=1,ncomp
             delta_eos = delta_eos + sedgex(i,j,k,comp)/rhobar(comp)
          end do

          do comp=1,ncomp
             sedgex(i,j,k,comp) = sedgex(i,j,k,comp) - w(comp)*(rhobar_sq/rhobar(comp))*delta_eos
          end do

          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains (rho,rho_1) instead of rho_i
             temp = sedgex(i,j,k,1)
             sedgex(i,j,k,1) = sedgex(i,j,k,1) + sedgex(i,j,k,2)
             sedgex(i,j,k,2) = temp
          end if
          
       enddo
       enddo
       enddo

       ! L2 projection: y-faces
       do k = ylo(3),yhi(3)
       do j = ylo(2),yhi(2) 
       do i = ylo(1),yhi(1) 
          
          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains rho_i instead of (rho,rho_1)
             temp = sedgey(i,j,k,1)
             sedgey(i,j,k,1) = sedgey(i,j,k,2)
             sedgey(i,j,k,2) = temp - sedgey(i,j,k,1)
          end if

          ! compute mass fractions, w_i = rho_i / rho
          temp = 0.d0
          do comp=1,ncomp
             temp = temp + sedgey(i,j,k,comp)
          end do
          do comp=1,ncomp
             w(comp) = sedgey(i,j,k,comp)/temp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do comp=1,ncomp
             rhobar_sq = rhobar_sq + w(comp)/rhobar(comp)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do comp=1,ncomp
             delta_eos = delta_eos + sedgey(i,j,k,comp)/rhobar(comp)
          end do

          do comp=1,ncomp
             sedgey(i,j,k,comp) = sedgey(i,j,k,comp) - w(comp)*(rhobar_sq/rhobar(comp))*delta_eos
          end do

          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains (rho,rho_1) instead of rho_i
             temp = sedgey(i,j,k,1)
             sedgey(i,j,k,1) = sedgey(i,j,k,1) + sedgey(i,j,k,2)
             sedgey(i,j,k,2) = temp
          end if
          
       enddo
       enddo
       enddo

       ! L2 projection: z-faces
       do k = zlo(3),zhi(3)
       do j = zlo(2),zhi(2) 
       do i = zlo(1),zhi(1) 
          
          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains rho_i instead of (rho,rho_1)
             temp = sedgez(i,j,k,1)
             sedgez(i,j,k,1) = sedgez(i,j,k,2)
             sedgez(i,j,k,2) = temp - sedgez(i,j,k,1)
          end if

          ! compute mass fractions, w_i = rho_i / rho
          temp = 0.d0
          do comp=1,ncomp
             temp = temp + sedgez(i,j,k,comp)
          end do
          do comp=1,ncomp
             w(comp) = sedgez(i,j,k,comp)/temp
          end do

          ! rhobar_sq = (sum_i (w_i/rhobar_i^2))^-1
          rhobar_sq = 0.d0
          do comp=1,ncomp
             rhobar_sq = rhobar_sq + w(comp)/rhobar(comp)**2
          end do
          rhobar_sq = 1.d0/rhobar_sq

          ! delta_eos = sum_l (rho_l/rhobar_l) - 1
          delta_eos = -1.d0
          do comp=1,ncomp
             delta_eos = delta_eos + sedgez(i,j,k,comp)/rhobar(comp)
          end do

          do comp=1,ncomp
             sedgez(i,j,k,comp) = sedgez(i,j,k,comp) - w(comp)*(rhobar_sq/rhobar(comp))*delta_eos
          end do

          if (proj_type .eq. 1) then
             ! overwrite sedge so it contains (rho,rho_1) instead of rho_i
             temp = sedgez(i,j,k,1)
             sedgez(i,j,k,1) = sedgez(i,j,k,1) + sedgez(i,j,k,2)
             sedgez(i,j,k,2) = temp
          end if
          
       enddo
       enddo
       enddo

    end if

  end subroutine bdsupdate_3d

  subroutine bdsfluxdiv_3d(lo,hi,s_update,ng_u,uadv,vadv,wadv,ng_v, &
                          sedgex,sedgey,sedgez,ng_e,dx,ncomp, &
                          tlo,thi)

    integer        ,intent(in   ) :: lo(:),hi(:),ng_u,ng_v,ng_e,ncomp
    integer        ,intent(in   ) :: tlo(:),thi(:)
    real(kind=dp_t),intent(inout) :: s_update(lo(1)-ng_u:,lo(2)-ng_u:,lo(3)-ng_u:,:)
    real(kind=dp_t),intent(in   ) ::     uadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     vadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(in   ) ::     wadv(lo(1)-ng_v:,lo(2)-ng_v:,lo(3)-ng_v:)
    real(kind=dp_t),intent(inout) ::   sedgex(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgey(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(inout) ::   sedgez(lo(1)-ng_e:,lo(2)-ng_e:,lo(3)-ng_e:,:)
    real(kind=dp_t),intent(in   ) :: dx(:)

    ! local variables
    integer comp,i,j,k

    real(kind=dp_t) :: w(ncomp), rhobar_sq, delta_eos, temp, dxinv

    ! advance solution
    ! conservative update

    dxinv=1.d0/dx(1) !assumed all dx equal
    

    do comp=1,ncomp
       do k = tlo(3),thi(3)
       do j = tlo(2),thi(2) 
       do i = tlo(1),thi(1) 
          s_update(i,j,k,comp) = s_update(i,j,k,comp) - (  &
               (sedgex(i+1,j,k,comp)*uadv(i+1,j,k)-sedgex(i,j,k,comp)*uadv(i,j,k))*dxinv +  &
               (sedgey(i,j+1,k,comp)*vadv(i,j+1,k)-sedgey(i,j,k,comp)*vadv(i,j,k))*dxinv + &
               (sedgez(i,j,k+1,comp)*wadv(i,j,k+1)-sedgez(i,j,k,comp)*wadv(i,j,k))*dxinv )
       enddo
       enddo
       enddo
    enddo

    end subroutine bdsfluxdiv_3d

  subroutine eval(s,slope,del,val)

    real(kind=dp_t), intent(in   ) :: s
    real(kind=dp_t), intent(in   ) :: slope(:)
    real(kind=dp_t), intent(in   ) :: del(:)
    real(kind=dp_t), intent(  out) :: val

    val = s + del(1)*slope(1) + del(2)*slope(2) + del(3)*slope(3) &
         + del(1)*del(2)*slope(4) + del(1)*del(3)*slope(5) + del(2)*del(3)*slope(6) &
         + del(1)*del(2)*del(3)*slope(7)

  end subroutine eval

  subroutine bds_quad(mla,umac,s,s_update,force,s_fc,dx,dt,scomp,ncomp, &
                      bc_comp,the_bc_tower,proj_type_in)
    ! modified for having the quadratic terms as well
    ! slxx and slyy are 2nd derivatives
    ! ave is the new constant for the polynomial

    type(ml_layout), intent(in   ) :: mla
    type(multifab) , intent(in   ) :: s(:)
    type(multifab) , intent(inout) :: s_update(:)
    type(multifab) , intent(in   ) :: force(:)
    type(multifab) , intent(in   ) :: s_fc(:,:)
    type(multifab) , intent(in   ) :: umac(:,:)
    real(kind=dp_t), intent(in   ) :: dx(:,:),dt
    integer        , intent(in   ) :: scomp, ncomp, bc_comp
    type(bc_tower) , intent(in   ) :: the_bc_tower
    integer        , intent(in   ), optional :: proj_type_in

    type(multifab) :: ave,slx,sly,slxy,slxx,slyy,sint,sc

    real(kind=dp_t), pointer :: uadvp(:,:,:,:)
    real(kind=dp_t), pointer :: vadvp(:,:,:,:)
    real(kind=dp_t), pointer ::   sop(:,:,:,:)
    real(kind=dp_t), pointer ::   sup(:,:,:,:)
    real(kind=dp_t), pointer ::    fp(:,:,:,:)
    real(kind=dp_t), pointer ::  avep(:,:,:,:)
    real(kind=dp_t), pointer ::  slxp(:,:,:,:)
    real(kind=dp_t), pointer ::  slyp(:,:,:,:)
    real(kind=dp_t), pointer :: slxyp(:,:,:,:)
    real(kind=dp_t), pointer :: slxxp(:,:,:,:)
    real(kind=dp_t), pointer :: slyyp(:,:,:,:)
    real(kind=dp_t), pointer ::   sip(:,:,:,:)
    real(kind=dp_t), pointer ::   scp(:,:,:,:)
    real(kind=dp_t), pointer :: spx(:,:,:,:)
    real(kind=dp_t), pointer :: spy(:,:,:,:)
    real(kind=dp_t), pointer :: sxp(:,:,:,:)
    real(kind=dp_t), pointer :: syp(:,:,:,:)

    integer :: dm,ng_s,ng_u,ng_o,ng_f,ng_e,ng_v,comp,bccomp,lev,i
    integer :: lo(2),hi(2),proj_type

    type(multifab) :: sedge(mla%dim)

    type(mfiter) :: mfi
    type(box) :: xnodalbox, ynodalbox
    integer :: xlo(mla%dim), xhi(mla%dim)
    integer :: ylo(mla%dim), yhi(mla%dim)
    type(box) :: tilebox
    integer :: tlo(mla%dim), thi(mla%dim)

    type(bl_prof_timer),save :: bpt

    call build(bpt,"bds_quad")

    ! L2 projection onto EOS?
    ! 0 = do nothing
    ! 1 = assume s holds rho and rho*c
    ! 2 = assume s holds all rho_i's
    if (present(proj_type_in)) then
       proj_type = proj_type_in
    else
       proj_type = 0
    end if

    ! Only worry about one level
    lev = 1
    dm = mla%dim

    do i=1,dm
       call multifab_build_edge(sedge(i),mla%la(1),ncomp,0,i)
    end do

    ! These all have one ghost cell and one component
    call multifab_build( ave,mla%la(lev),1,1)
    call multifab_build( slx,mla%la(lev),1,1)
    call multifab_build( sly,mla%la(lev),1,1)
    call multifab_build(slxy,mla%la(lev),1,1)
    call multifab_build(slxx,mla%la(lev),1,1)
    call multifab_build(slyy,mla%la(lev),1,1)

    ! These has two ghost cells and one component
    call multifab_build(sint,mla%la(lev),1,2)

    ! This has one ghost cell and four components
    call multifab_build(  sc,mla%la(lev),4,1)

    ng_s = s(1)%ng 
    ng_u = s_update(1)%ng
    ng_o = force(1)%ng
    ng_f = s_fc(1,1)%ng
    ng_e = sedge(1)%ng
    ng_v = umac(1,1)%ng

    do i = 1, nfabs(s(lev))
       uadvp => dataptr(umac(lev,1), i)
       vadvp => dataptr(umac(lev,2), i)
       sop   => dataptr(s(lev) , i)
       sup   => dataptr(s_update(lev), i)
       fp    => dataptr(force(lev), i)
       spx => dataptr(s_fc(lev,1), i)
       spy => dataptr(s_fc(lev,2), i)
       avep  => dataptr(ave , i)
       slxp  => dataptr(slx , i)
       slyp  => dataptr(sly , i)
       slxyp => dataptr(slxy, i)
       slxxp  => dataptr(slxx , i)
       slyyp  => dataptr(slyy , i)
       sip   => dataptr(sint, i)
       scp   => dataptr(sc  , i)
       sxp => dataptr(sedge(1), i)
       syp => dataptr(sedge(2), i)
       lo =  lwb(get_box(s(lev), i))
       hi =  upb(get_box(s(lev), i))
       do comp = scomp, scomp+ncomp-1
          bccomp = bc_comp+comp-scomp
          select case (dm)
          case (2)

             if (advection_type .eq. 3) then

                call bdsslope_quad_2d(lo, hi, sop(:,:,1,comp), ng_s, &
                                      avep(:,:,1,1), slxp(:,:,1,1), slyp(:,:,1,1), &
                                      slxyp(:,:,1,1), slxxp(:,:,1,1), slyyp(:,:,1,1), &
                                      sip(:,:,1,1), scp(:,:,1,:), dx(lev,:)) 

             else if (advection_type .eq. 4) then

                call bdsslope_quad_lim_2d(lo, hi, sop(:,:,1,comp), ng_s, &
                                          avep(:,:,1,1), slxp(:,:,1,1), slyp(:,:,1,1), &
                                          slxyp(:,:,1,1), slxxp(:,:,1,1), slyyp(:,:,1,1), &
                                          sip(:,:,1,1), scp(:,:,1,:), dx(lev,:)) 

             end if

             call  bdsconc_quad_2d(lo, hi, sup(:,:,1,comp), ng_u, &
                                   fp(:,:,1,comp), ng_o, &
                                   avep(:,:,1,1), slxp(:,:,1,1), slyp(:,:,1,1), slxyp(:,:,1,1), &
                                   slxxp(:,:,1,1), slyyp(:,:,1,1), &
                                   uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, dx(lev,:), dt, &
                                   spx(:,:,1,comp), spy(:,:,1,comp), ng_f, &
                                   sxp(:,:,1,comp-scomp+1), &
                                   syp(:,:,1,comp-scomp+1), ng_e, &
                                   the_bc_tower%bc_tower_array(1)%adv_bc_level_array(i,:,:,bccomp))
          case (3)
             call parallel_abort("quadratic BDS advection not supported in 3D")  
          end select
       end do
    end do

    ! do L2 projection and increment s_update
    
    !$omp parallel private(mfi,i,xnodalbox,xlo,xhi,ynodalbox,ylo,yhi) &
    !$omp private(sup,sxp,syp,uadvp,vadvp,lo,hi)
    do i = 1, nfabs(s(lev))
       sup => dataptr(s_update(lev), i)
       sxp => dataptr(sedge(1), i)
       syp => dataptr(sedge(2), i)
       uadvp  => dataptr(umac(lev,1), i)
       vadvp  => dataptr(umac(lev,2), i)
       lo =  lwb(get_box(s(lev), i))
       hi =  upb(get_box(s(lev), i))
       select case (dm)
       case (2)
          call bdsupdate_2d(lo, hi, &
                            uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, &
                            sxp(:,:,1,:), syp(:,:,1,:), ng_e, &
                            dx(lev,:), ncomp, proj_type,xlo,xhi,ylo,yhi)
       case (3)
          call parallel_abort("quadratic BDS advection not supported in 3D")  
       end select
    end do ! end loop over fabs
    !$omp end parallel

    ! do L2 projection and increment s_update

    !$omp parallel private(mfi,i,tilebox,tlo,thi) &
    !$omp private(sup,sxp,syp) &
    !$omp private(uadvp,vadvp,lo,hi)
    call mfiter_build(mfi, s_update(lev), tiling=.true.)
    do while (more_tile(mfi))
       i = get_fab_index(mfi)

       tilebox = get_tilebox(mfi)
       tlo = lwb(tilebox)
       thi = upb(tilebox)

       sup => dataptr(s_update(lev), i)
       sxp => dataptr(sedge(1), i)
       syp => dataptr(sedge(2), i)
       uadvp  => dataptr(umac(lev,1), i)
       vadvp  => dataptr(umac(lev,2), i)
       lo =  lwb(get_box(s(lev), i))
       hi =  upb(get_box(s(lev), i))
       select case (dm)
       case (2)
          call bdsfluxdiv_2d(lo, hi, &
                            sup(:,:,1,scomp:), ng_u, &
                            uadvp(:,:,1,1), vadvp(:,:,1,1), ng_v, &
                            sxp(:,:,1,:), syp(:,:,1,:), ng_e, &
                            dx(lev,:), ncomp,tlo,thi)
       case (3)
          call parallel_abort("quadratic BDS advection not supported in 3D")  
       end select
    end do ! end loop over fabs
    !$omp end parallel

    call multifab_destroy(ave)
    call multifab_destroy(slx)
    call multifab_destroy(sly)
    call multifab_destroy(slxy)
    call multifab_destroy(slxx)
    call multifab_destroy(slyy)
    call multifab_destroy(sint)
    call multifab_destroy(sc)
    do i=1,dm
       call multifab_destroy(sedge(i))
    end do

    call destroy(bpt)

  end subroutine bds_quad

  subroutine bdsslope_quad_2d(lo,hi,s,ng_s,ave,slx,sly,slxy,slxx,slyy,sint,sc,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng_s
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng_s:,lo(2)-ng_s:)
    real(kind=dp_t), intent(inout) ::  ave(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) ::  slx(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) ::  sly(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slxy(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slxx(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slyy(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: sint(lo(1)-2:,lo(2)-2:)
    real(kind=dp_t), intent(inout) ::   sc(lo(1)-1:,lo(2)-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)

    real(kind=dp_t), allocatable :: diff(:,:)
    real(kind=dp_t), allocatable :: smin(:,:)
    real(kind=dp_t), allocatable :: smax(:,:)
    real(kind=dp_t), allocatable :: sumdif(:)
    real(kind=dp_t), allocatable :: sgndif(:)
    integer        , allocatable :: kdp(:)

    real(kind=dp_t) :: hx,hy,sumloc,redfac,redmax,div,hxinv,hyinv
    real(kind=dp_t) :: eps
    integer         :: inc1, inc2, inc3, inc4
    integer         :: i,j,k,ll,is,ie,js,je

    allocate(  diff(lo(1)-1:hi(1)+1,4))
    allocate(  smin(lo(1)-1:hi(1)+1,4))
    allocate(  smax(lo(1)-1:hi(1)+1,4))
    allocate(sumdif(lo(1)-1:hi(1)+1  ))
    allocate(sgndif(lo(1)-1:hi(1)+1  ))
    allocate(   kdp(lo(1)-1:hi(1)+1  ))

    hx = dx(1)
    hy = dx(2)
    hxinv = 1/hx
    hyinv = 1/hy
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    eps = 1.d-10

    ! estimate corner values and calculate slopes out of it

    do i = is-2,ie+1
       do j = js-2,je+1
          sint(i,j) = (s(i-1,j-1) + s(i-1,j+2) + s(i+2,j-1) + s(i+2,j+2) &
               - 7.d0*(s(i-1,j  ) + s(i-1,j+1) + s(i  ,j-1) + s(i+1,j-1) + & 
                       s(i  ,j+2) + s(i+1,j+2) + s(i+2,j  ) + s(i+2,j+1)) +  &
                49.d0*(s(i  ,j  ) + s(i+1,j  ) + s(i  ,j+1) + s(i+1,j+1)) ) / 144.d0
       enddo
    enddo

    do j = js-1,je+1
       do i = is-1,ie+1 

          slx(i,j) = 0.5d0*(sint(i  ,j) + sint(i  ,j-1) - &
                            sint(i-1,j) - sint(i-1,j-1) ) * hxinv
          sly(i,j) = 0.5d0*(sint(i  ,j) - sint(i  ,j-1) + &
                            sint(i-1,j) - sint(i-1,j-1) ) * hyinv
          slxy(i,j) = (sint(i,j  ) - sint(i  ,j-1) - &
                       sint(i-1,j) + sint(i-1,j-1) ) * (hxinv*hyinv)

       enddo
    enddo

    ! estimate 2nd derivatives at cell center (slxx and slyy are 1/2 of that) 
    ! and adjust constant to fit cell average

    do j = js-1,je+1
       do i = is-1,ie+1 

          slxx(i,j) = 0.5d0*( - s(i-2,j) + 12.d0*s(i-1,j) - 22.d0*s(i,j)  &
               + 12.d0*s(i+1,j) - s(i+2,j) ) / (8.d0*hx**2)
          slyy(i,j) = 0.5d0*( - s(i,j-2) + 12.d0*s(i,j-1) - 22.d0*s(i,j)  &
               + 12.d0*s(i,j+1) - s(i,j+2) ) / (8.d0*hy**2)
            
          ave(i,j) = s(i,j) - ( slxx(i,j)*hx**2 + slyy(i,j)*hy**2 ) / 12.d0

       enddo
    enddo
      
    deallocate(diff,smin,smax,sumdif,sgndif,kdp)

  end subroutine bdsslope_quad_2d

  subroutine bdsslope_quad_lim_2d(lo,hi,s,ng,ave,slx,sly,slxy,slxx,slyy,sint,sc,dx)

    integer        , intent(in   ) :: lo(:), hi(:), ng
    real(kind=dp_t), intent(in   ) ::    s(lo(1)-ng:,lo(2)-ng:)
    real(kind=dp_t), intent(inout) ::  ave(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) ::  slx(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) ::  sly(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slxy(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slxx(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: slyy(lo(1)-1:,lo(2)-1:)
    real(kind=dp_t), intent(inout) :: sint(lo(1)-2:,lo(2)-2:)
    real(kind=dp_t), intent(inout) ::   sc(lo(1)-1:,lo(2)-1:,:)
    real(kind=dp_t), intent(in   ) :: dx(:)


    real(kind=dp_t), allocatable :: diff(:,:)
    real(kind=dp_t), allocatable :: smin(:,:)
    real(kind=dp_t), allocatable :: smax(:,:)
    real(kind=dp_t), allocatable :: sumdif(:)
    real(kind=dp_t), allocatable :: sgndif(:)
    integer        , allocatable :: kdp(:)
    logical        , allocatable :: test_free(:,:)
    real(kind=dp_t), allocatable :: adjust_goal(:,:)
    real(kind=dp_t), allocatable :: slxx_extr(:,:)
    real(kind=dp_t), allocatable :: slyy_extr(:,:)


    real(kind=dp_t) :: hx,hy,sumloc,redfac,redmax,div,cmp,cmp_x,cmp_y,my_sn,min_x,min_y
    real(kind=dp_t) :: eps, pos, my_val,my_min,my_max,hxinv,hyinv
    integer         :: inc1, inc2, inc3, inc4, ind_x1, ind_y1,ind_x2, ind_y2, limit_x, limit_y 
    integer         :: edge_xh, edge_xl, edge_yh, edge_yl, test_corners
    integer         :: i,j,k,ll,is,ie,js,je

    allocate(  diff(lo(1)-1:hi(1)+1,4))
    allocate(  smin(lo(1)-1:hi(1)+1,4))
    allocate(  smax(lo(1)-1:hi(1)+1,4))
    allocate(sumdif(lo(1)-1:hi(1)+1  ))
    allocate(sgndif(lo(1)-1:hi(1)+1  ))
    allocate(   kdp(lo(1)-1:hi(1)+1  ))
    allocate(  test_free(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(  adjust_goal(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1))
    allocate(  slxx_extr(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2))
    allocate(  slyy_extr(lo(1)-2:hi(1)+2,lo(2)-2:hi(2)+2))


    hx = dx(1)
    hy = dx(2)
    hxinv = 1/hx
    hyinv = 1/hy
    is = lo(1)
    ie = hi(1)
    js = lo(2)
    je = hi(2)

    eps = 1.d-11

    ! extrapolate corner values and second deriavtives

    do i = is-2,ie+1
       do j = js-2,je+1
          sint(i,j) = (s(i-1,j-1) + s(i-1,j+2) + s(i+2,j-1) + s(i+2,j+2) &
               - 7.d0*(s(i-1,j  ) + s(i-1,j+1) + s(i  ,j-1) + s(i+1,j-1) + & 
               s(i  ,j+2) + s(i+1,j+2) + s(i+2,j  ) + s(i+2,j+1)) +  &
               49.d0*(s(i  ,j  ) + s(i+1,j  ) + s(i  ,j+1) + s(i+1,j+1)) ) / 144.d0
       end do
    end do

    do j = js-2,je+2
       do i = is-2,ie+2 
          slxx_extr(i,j) = 0.5d0*( - s(i-2,j) + 12.d0*s(i-1,j) - 22.d0*s(i,j)  &
               + 12.d0*s(i+1,j) - s(i+2,j) ) / (8.d0*hx**2)
          slyy_extr(i,j) = 0.5d0*( - s(i,j-2) + 12.d0*s(i,j-1) - 22.d0*s(i,j)  &
               + 12.d0*s(i,j+1) - s(i,j+2) ) / (8.d0*hy**2)
       enddo
    enddo


    ! if all extrapolated values are on one side of s(i,j), construct a constant on that cell
    ! set test_free to false so that the polynomial is not changed in what follows
    do j = js-1,je+1
       do i = is-1,ie+1
          if (  min((sint(i,j)- s(i,j))*(sint(i-1,j)-s(i,j)) , &
               (sint(i,j)- s(i,j))*(sint(i,j-1)-s(i,j)), &
               (sint(i,j)- s(i,j))*(sint(i-1,j-1)-s(i,j)) ) .ge. 0.d0  ) then
             test_free(i,j) = .false.
             slx(i,j) = 0.d0
             sly(i,j) = 0.d0
             slxy(i,j) = 0.d0
             slxx(i,j) = 0.d0
             slyy(i,j) = 0.d0
             ave(i,j) = s(i,j)
          else
             test_free(i,j) = .true.
          end if
       end do
    end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! first round: try if unlimited polynomial passes the test
    ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    do j = js-1,je+1
       do i = is-1,ie+1 

          if (test_free(i,j) ) then

             ! calculate slopes out of sint

             slx(i,j) = 0.5d0*(sint(i  ,j) + sint(i  ,j-1) - &
                  sint(i-1,j) - sint(i-1,j-1) ) * hxinv
             sly(i,j) = 0.5d0*(sint(i  ,j) - sint(i  ,j-1) + &
                  sint(i-1,j) - sint(i-1,j-1) ) * hyinv
             slxy(i,j) = (sint(i,j  ) - sint(i  ,j-1) - &
                  sint(i-1,j) + sint(i-1,j-1) ) * (hxinv*hyinv)

             ! limit slxx, slyy such that no true interior min/max

             slxx(i,j) = slxx_extr(i,j)

             slyy(i,j) = slyy_extr(i,j)

             ! make sure that there is no interior min/max
             ind_x1=0;
             ind_x2=0;
             ind_y1=0;
             ind_y2=0;
             if (sign(1.d0,-slxy(i,j)*hy/2.-slx(i,j)) * sign(1.d0,slxy(i,j)*hy/2.-slx(i,j)) .lt. 0.d0 ) then
                ind_x1 = 1;
             else
                cmp_x = min(abs(slxy(i,j)*hy/2.+slx(i,j)) , abs(slxy(i,j)*hy/2.-slx(i,j)) )
                if (cmp_x .lt. abs(slxx(i,j))*hx) then ! ie would create interior extrema
                   ind_x2 = 1;
                end if
             end if

             if (sign(1.d0,-slxy(i,j)*hx/2.-sly(i,j)) * sign(1.d0,slxy(i,j)*hx/2.-sly(i,j)) .lt. 0.d0 ) then
                ind_y1 = 1;
             else
                cmp_y = min(abs(slxy(i,j)*hx/2.+sly(i,j)) , abs(slxy(i,j)*hx/2.-sly(i,j)) )
                if (cmp_y .lt. abs(slyy(i,j))*hy) then ! ie would create interior extrema
                   ind_y2 = 1;
                end if
             end if

             ! set slxx, slyy=0 if real interior min/max (ie both x-der and y-der are 0)
             if (ind_x1 .eq. 1 .and. ind_y1 .eq. 1) then
                slxx(i,j) = 0.d0;
                slyy(i,j) = 0.d0;
                !             print *,'set slxx,slyy=0,i,j',i,j
             end if
             if (ind_x1 .eq. 1 .and. ind_y2 .eq. 1) then
                slxx(i,j) = 0.d0;
                slyy(i,j) = sign(1.0d0,slyy_extr(i,j)) * cmp_y * hyinv;
             end if
             if (ind_x2 .eq. 1 .and. ind_y1 .eq. 1) then
                slxx(i,j) =  sign(1.0d0,slxx_extr(i,j)) * cmp_x * hxinv;
                slyy(i,j) = 0.d0;
             end if
             if (ind_x2 .eq. 1 .and. ind_y2 .eq. 1) then
                slxx(i,j) =  sign(1.0d0,slxx_extr(i,j)) * cmp_x * hxinv;
                slyy(i,j) =  sign(1.0d0,slyy_extr(i,j)) * cmp_y * hyinv;
             end if


             ! Then adjust constant
             ave(i,j) = s(i,j) - ( slxx(i,j)*hx**2 + slyy(i,j)*hy**2 ) / 12.d0

!!!!!!!!!!!!!!!!!!!!!!!!
             ! still have to check for extrema on the boundary
!!!!!!!!!!!!!!!!!!!!!!!!

             smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

             sc(i,j,4) = s(i,j) + 0.5d0*(hx*slx(i,j) + hy*sly(i,j))  &
                  + 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,3) = s(i,j) + 0.5d0*(hx*slx(i,j) - hy*sly(i,j))  &
                  - 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,2) = s(i,j) - 0.5d0*(hx*slx(i,j) - hy*sly(i,j)) &
                  - 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,1) = s(i,j) - 0.5d0*(hx*slx(i,j) + hy*sly(i,j)) &
                  + 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.

             if (sc(i,j,4) .le. smax(i,4) .and. sc(i,j,4) .ge. smin(i,4)  .and. &
                  sc(i,j,3) .le. smax(i,3) .and. sc(i,j,3) .ge. smin(i,3) .and. &
                  sc(i,j,2) .le. smax(i,2) .and. sc(i,j,2) .ge. smin(i,2)  .and. &
                  sc(i,j,1) .le. smax(i,1) .and. sc(i,j,1) .ge. smin(i,1)) then
                ! keep that polynomial (if it passes the edge tests as well)
                test_corners = 1
             else
                test_corners = 0;
                if (j .ge. js .and. j .le. je .and. i .ge. is .and. i .le. ie) then
                   !                print *,'corners dont pass test,i,j',i,j
                end if
             end if

             ! upper y-edge
             if (abs(slxy(i,j)*hy/2.+slx(i,j)) .lt. abs(hx*slxx(i,j)) ) then
                ! have interior min/max
                pos = - (slxy(i,j)*hy/2.+slx(i,j))/(2*slxx(i,j))
                my_val = ave(i,j) + pos*slx(i,j) + 0.5d0*hy*sly(i,j)  &
                     + 0.5d0*pos*hy*slxy(i,j) &
                     + pos**2*slxx(i,j) + hy**2*slyy(i,j)/4.
                if (pos .ge. 0.d0) then
                   my_min = min(s(i,j),s(i,j+1),s(i+1,j),s(i+1,j+1))
                   my_max = max(s(i,j),s(i,j+1),s(i+1,j),s(i+1,j+1))
                else
                   my_min = min(s(i,j),s(i,j+1),s(i-1,j),s(i-1,j+1))
                   my_max = max(s(i,j),s(i,j+1),s(i-1,j),s(i-1,j+1))
                end if
                if (my_val .ge. my_min .and. my_val .le. my_max) then
                   edge_yh = 1;
                   limit_x = 0;
                else
                   edge_yh = 0;
                   limit_x = 1;
                   if (j .ge. js .and. j .le. je .and. i .ge. is .and. i .le. ie) then
                      !                   print *,'edge_yh doesnt pass test,i,j',i,j
                   end if
                end if
             else 
                edge_yh = 1;
             end if

             ! lower y-edge
             if (abs(-slxy(i,j)*hy/2.+slx(i,j)) .lt. abs(hx*slxx(i,j)) ) then
                ! have interior min/max
                pos = - (-slxy(i,j)*hy/2.+slx(i,j))/(2*slxx(i,j))
                my_val = ave(i,j) + pos*slx(i,j) - 0.5d0*hy*sly(i,j)  &
                     - 0.5d0*pos*hy*slxy(i,j) &
                     + pos**2*slxx(i,j) + hy**2*slyy(i,j)/4.
                if (pos .gt. 0.d0) then
                   my_min = min(s(i,j),s(i,j-1),s(i+1,j),s(i+1,j-1))
                   my_max = max(s(i,j),s(i,j-1),s(i+1,j),s(i+1,j-1))
                else
                   my_min = min(s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1))
                   my_max = max(s(i,j),s(i,j-1),s(i-1,j),s(i-1,j-1))
                end if
                if (my_val .ge. my_min .and. my_val .le. my_max) then
                   edge_yl = 1;
                else
                   edge_yl = 0;
                   limit_x = 1;
                   if (j .ge. js .and. j .le. je .and. i .ge. is .and. i .le. ie) then
                      !                   print *,'edge_yl doesnt pass test,i,j',i,j
                   end if
                end if
             else 
                edge_yl = 1;
             end if

             ! upper x-edge
             if (abs(slxy(i,j)*hx/2.+sly(i,j)) .lt. abs(hy*slyy(i,j)) ) then
                ! have interior min/max
                pos = - (slxy(i,j)*hx/2.+sly(i,j))/(2*slyy(i,j))
                my_val = ave(i,j) + 0.5*hx*slx(i,j) + pos*sly(i,j)  &
                     + 0.5d0*hx*pos*slxy(i,j) &
                     + hx**2*slxx(i,j)/4. + pos**2*slyy(i,j)
                if (pos .ge. 0.d0) then
                   my_min = min(s(i,j),s(i+1,j),s(i,j+1),s(i+1,j+1))
                   my_max = max(s(i,j),s(i+1,j),s(i,j+1),s(i+1,j+1))
                else
                   my_min = min(s(i,j),s(i+1,j),s(i,j-1),s(i+1,j-1))
                   my_max = max(s(i,j),s(i+1,j),s(i,j-1),s(i+1,j-1))
                end if
                if (my_val .ge. my_min .and. my_val .le. my_max) then
                   edge_xh = 1;
                   limit_y = 0;
                else
                   edge_xh = 0;
                   limit_y = 1;
                   if (j .ge. js .and. j .le. je .and. i .ge. is .and. i .le. ie) then
                      !                   print *,'edge_xh doesnt pass test,i,j',i,j
                   end if
                end if
             else 
                edge_xh = 1;
             end if

             ! lower x-edge
             if (abs(-slxy(i,j)*hx/2.+sly(i,j)) .lt. abs(hy*slyy(i,j)) ) then
                ! have interior min/max
                pos = - (-slxy(i,j)*hx/2.+sly(i,j))/(2*slyy(i,j))
                my_val = ave(i,j) - 0.5*hx*slx(i,j) + pos*sly(i,j)  &
                     - 0.5d0*hx*pos*slxy(i,j) &
                     + hx**2*slxx(i,j)/4. + pos**2*slyy(i,j)
                if (pos .ge. 0.d0) then
                   my_min = min(s(i,j),s(i-1,j),s(i,j+1),s(i-1,j+1))
                   my_max = max(s(i,j),s(i-1,j),s(i,j+1),s(i-1,j+1))
                else
                   my_min = min(s(i,j),s(i-1,j),s(i,j-1),s(i-1,j-1))
                   my_max = max(s(i,j),s(i-1,j),s(i,j-1),s(i-1,j-1))
                end if
                if (my_val .ge. my_min .and. my_val .le. my_max) then
                   edge_xl = 1;
                else
                   edge_xl = 0;
                   limit_y = 1;
                   if (j .ge. js .and. j .le. je .and. i .ge. is .and. i .le. ie) then
                      !                   print *,'edge_xl doesnt pass test,i,j',i,j
                   end if
                end if
             else 
                edge_xl = 1;
             end if

             if (test_corners .eq. 1 .and. edge_yh .eq. 1 .and. edge_yl .eq. 1 &
                  .and. edge_xh .eq. 1 .and. edge_xl .eq. 1) then
                ! keep that polynomial
                test_free(i,j) = .false.
                !             print *,'success in first round,i,j',i,j

             else ! limit slxx and/or slyy and check again

!!!!!!!! begin else

                if (sign(1.d0,-slxy(i,j)*hy/2.-slx(i,j)) * sign(1.d0,slxy(i,j)*hy/2.-slx(i,j)) .lt. 0.d0 ) then
                   slxx(i,j) = 0.d0
                else
                   cmp = min(abs(slxy(i,j)*hy/2.+slx(i,j)) , abs(slxy(i,j)*hy/2.-slx(i,j)) )
                   if (cmp .lt. abs(slxx(i,j))*hx) then ! ie would create interior extrema
                      slxx(i,j) = sign(1.0d0,slxx_extr(i,j)) * cmp * hxinv
                   end if
                end if

                if (sign(1.d0,-slxy(i,j)*hx/2.-sly(i,j)) * sign(1.d0,slxy(i,j)*hx/2.-sly(i,j)) .lt. 0.d0 ) then
                   slyy(i,j) = 0.d0
                else
                   cmp = min(abs(slxy(i,j)*hx/2.+sly(i,j)) , abs(slxy(i,j)*hx/2.-sly(i,j)) )
                   if (cmp .lt. abs(slyy(i,j))*hy) then ! ie would create interior extrema
                      slyy(i,j) = sign(1.0d0,slyy_extr(i,j)) * cmp * hyinv
                   end if
                end if

                ave(i,j) = s(i,j) - ( slxx(i,j)*hx**2 + slyy(i,j)*hy**2 ) / 12.d0

                ! Test again

                sc(i,j,4) = s(i,j) + 0.5d0*(hx*slx(i,j) + hy*sly(i,j))  &
                     + 0.25d0*hx*hy*slxy(i,j) &
                     + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
                sc(i,j,3) = s(i,j) + 0.5d0*(hx*slx(i,j) - hy*sly(i,j))  &
                     - 0.25d0*hx*hy*slxy(i,j) &
                     + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
                sc(i,j,2) = s(i,j) - 0.5d0*(hx*slx(i,j) - hy*sly(i,j)) &
                     - 0.25d0*hx*hy*slxy(i,j) &
                     + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
                sc(i,j,1) = s(i,j) - 0.5d0*(hx*slx(i,j) + hy*sly(i,j)) &
                     + 0.25d0*hx*hy*slxy(i,j) &
                     + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.

                if (sc(i,j,4) .le. smax(i,4) .and. sc(i,j,4) .ge. smin(i,4)  .and. &
                     sc(i,j,3) .le. smax(i,3) .and. sc(i,j,3) .ge. smin(i,3) .and. &
                     sc(i,j,2) .le. smax(i,2) .and. sc(i,j,2) .ge. smin(i,2)  .and. &
                     sc(i,j,1) .le. smax(i,1) .and. sc(i,j,1) .ge. smin(i,1)) then
                   ! keep that polynomial 
                   test_free(i,j) = .false.
                end if
!!!!! end else 

             end if

          end if

       enddo
    enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! if unlimited slopes not successful: go back to routine bds_quadr_monotonelimiting
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! NOTE: in the following we only work on cells for which test_free .eq. .true.

    ! we already extrapolated corner values
    ! and calculated slopes


    do j = js-1,je+1

       ! make sure we are in min/max bounds
       do i = is-1,ie+1
          if (test_free(i,j) ) then
             smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

             sc(i,j,4) = s(i,j) + 0.5d0*(hx*slx(i,j) + hy*sly(i,j))  &
                  + 0.25d0*hx*hy*slxy(i,j)
             sc(i,j,3) = s(i,j) + 0.5d0*(hx*slx(i,j) - hy*sly(i,j))  &
                  - 0.25d0*hx*hy*slxy(i,j)
             sc(i,j,2) = s(i,j) - 0.5d0*(hx*slx(i,j) - hy*sly(i,j)) &
                  - 0.25d0*hx*hy*slxy(i,j)
             sc(i,j,1) = s(i,j) - 0.5d0*(hx*slx(i,j) + hy*sly(i,j)) &
                  + 0.25d0*hx*hy*slxy(i,j)

             sc(i,j,4) = max(min(sc(i,j,4), smax(i,4)), smin(i,4))
             sc(i,j,3) = max(min(sc(i,j,3), smax(i,3)), smin(i,3))
             sc(i,j,2) = max(min(sc(i,j,2), smax(i,2)), smin(i,2))
             sc(i,j,1) = max(min(sc(i,j,1), smax(i,1)), smin(i,1))
          end if
       end do

       ! make average of corner values fit the cell average 
       do ll = 1,3 
          do i = is-1,ie+1 
             if (test_free(i,j) ) then
                sumloc = 0.25d0*(sc(i,j,4) + sc(i,j,3) +  &
                     sc(i,j,2) + sc(i,j,1))
                sumdif(i) = (sumloc - s(i,j))*4.d0
                sgndif(i) = sign(1.d0,sumdif(i))

                diff(i,4) = (sc(i,j,4) - s(i,j))*sgndif(i)
                diff(i,3) = (sc(i,j,3) - s(i,j))*sgndif(i)
                diff(i,2) = (sc(i,j,2) - s(i,j))*sgndif(i)
                diff(i,1) = (sc(i,j,1) - s(i,j))*sgndif(i)

                if (diff(i,1) .gt. eps) then
                   inc1 = 1
                else
                   inc1 = 0
                end if

                if (diff(i,2) .gt. eps) then
                   inc2 = 1
                else
                   inc2 = 0
                end if

                if (diff(i,3) .gt. eps) then
                   inc3 = 1
                else
                   inc3 = 0
                end if

                if (diff(i,4) .gt. eps) then
                   inc4 = 1
                else
                   inc4 = 0
                end if

                kdp(i) = inc1 + inc2 + inc3 + inc4
             end if
          enddo

          do k = 1,4 
             do i = is-1,ie+1 
                if (test_free(i,j) ) then
                   if (kdp(i).lt.1) then 
                      div = 1.d0
                   else
                      div = dble(kdp(i))
                   end if

                   if (diff(i,k).gt.eps) then
                      redfac = sumdif(i)*sgndif(i)/div
                      kdp(i) = kdp(i)-1
                   else
                      redfac = 0.d0
                   end if

                   if (sgndif(i) .gt. 0.d0) then
                      redmax = sc(i,j,k) - smin(i,k)
                   else
                      redmax = smax(i,k) - sc(i,j,k)
                   end if

                   redfac = min(redfac,redmax)
                   sumdif(i) = sumdif(i) - redfac*sgndif(i)
                   sc(i,j,k) = sc(i,j,k) - redfac*sgndif(i)
                end if
             enddo
          enddo

       enddo

       ! recalculate slopes
       do i = is-1,ie+1 
          if (test_free(i,j) ) then
             slx(i,j) = 0.5d0*(sc(i,j,4) + sc(i,j,3) -  &
                  sc(i,j,1) - sc(i,j,2))/hx
             sly(i,j) = 0.5d0*(sc(i,j,4) + sc(i,j,2) -  &
                  sc(i,j,1) - sc(i,j,3))/hy
             slxy(i,j) = ( sc(i,j,1) + sc(i,j,4) &
                  -sc(i,j,2) - sc(i,j,3) ) * (hxinv*hyinv)

          end if
       enddo

    end do  ! end of j-loop


    ! now add quadratic terms

    ! limit the second derivatives

    do j = js-1,je+1
       do i = is-1,ie+1 

          if (test_free(i,j)) then
             slxx(i,j) = slxx_extr(i,j)   

             ! make sure that reconstructed function is monotone

             if (sign(1.d0,-slxy(i,j)*hy/2.-slx(i,j)) * sign(1.d0,slxy(i,j)*hy/2.-slx(i,j)) .lt. 0.d0 ) then
                slxx(i,j) = 0.d0
             else
                cmp = min(abs(slxy(i,j)*hy/2.+slx(i,j)) , abs(slxy(i,j)*hy/2.-slx(i,j)) )
                if (cmp .lt. abs(slxx(i,j))*hx) then ! ie would create interior extrema
                   slxx(i,j) = sign(1.0d0,slxx_extr(i,j)) * cmp * hxinv
                end if
             end if

             ! do the same for slyy

             slyy(i,j) = slyy_extr(i,j)

             if (sign(1.d0,-slxy(i,j)*hx/2.-sly(i,j)) * sign(1.d0,slxy(i,j)*hx/2.-sly(i,j)) .lt. 0.d0 ) then
                slyy(i,j) = 0.d0
             else
                cmp = min(abs(slxy(i,j)*hx/2.+sly(i,j)) , abs(slxy(i,j)*hx/2.-sly(i,j)) )
                if (cmp .lt. abs(slyy(i,j))*hy) then ! ie would create interior extrema
                   slyy(i,j) = sign(1.0d0,slyy_extr(i,j)) * cmp * hyinv
                end if
             end if

             ! Then adjust constant

             ave(i,j) = s(i,j) - ( slxx(i,j)*hx**2 + slyy(i,j)*hy**2 ) / 12.d0

             ! Check whether we are still in max/min bounds:
             ! if not set slxx, slyy to 0
             smin(i,4) = min(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smax(i,4) = max(s(i,j), s(i+1,j), s(i,j+1), s(i+1,j+1))
             smin(i,3) = min(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smax(i,3) = max(s(i,j), s(i+1,j), s(i,j-1), s(i+1,j-1))
             smin(i,2) = min(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smax(i,2) = max(s(i,j), s(i-1,j), s(i,j+1), s(i-1,j+1))
             smin(i,1) = min(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))
             smax(i,1) = max(s(i,j), s(i-1,j), s(i,j-1), s(i-1,j-1))

             sc(i,j,4) = s(i,j) + 0.5d0*(hx*slx(i,j) + hy*sly(i,j))  &
                  + 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,3) = s(i,j) + 0.5d0*(hx*slx(i,j) - hy*sly(i,j))  &
                  - 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,2) = s(i,j) - 0.5d0*(hx*slx(i,j) - hy*sly(i,j)) &
                  - 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.
             sc(i,j,1) = s(i,j) - 0.5d0*(hx*slx(i,j) + hy*sly(i,j)) &
                  + 0.25d0*hx*hy*slxy(i,j) &
                  + hx**2*slxx(i,j)/6. + hy**2*slyy(i,j)/6.

             if (sc(i,j,4) .gt. smax(i,4)+1.e-10 .or. sc(i,j,4) .lt. smin(i,4)-1.e-10) then
                slxx(i,j) = 0
                slyy(i,j) = 0
                ave(i,j) = s(i,j)
                !            print *,'dropped order to ensure max/min,i,j',i,j
             end if
             if (sc(i,j,3) .gt. smax(i,3)+1.e-10 .or. sc(i,j,3) .lt. smin(i,3)-1.e-10) then
                slxx(i,j) = 0
                slyy(i,j) = 0
                ave(i,j) = s(i,j)
                !            print *,'dropped order to ensure max/min,i,j',i,j
             end if
             if (sc(i,j,2) .gt. smax(i,2)+1.e-10 .or. sc(i,j,2) .lt. smin(i,2)-1.e-10) then
                slxx(i,j) = 0
                slyy(i,j) = 0
                ave(i,j) = s(i,j)
                !            print *,'dropped order to ensure max/min,i,j',i,j
             end if
             if (sc(i,j,1) .gt. smax(i,1)+1.e-10 .or. sc(i,j,1) .lt. smin(i,1)-1.e-10) then
                slxx(i,j) = 0
                slyy(i,j) = 0
                ave(i,j) = s(i,j)
                !            print *,'dropped order to ensure max/min,i,j',i,j
             end if

          end if

       enddo
    enddo

    deallocate(diff,smin,smax,sumdif,sgndif,kdp,test_free,adjust_goal)

  end subroutine bdsslope_quad_lim_2d

  ! ***********************************************
  ! start of routine bdsconc_2d
  ! ***********************************************

  subroutine bdsconc_quad_2d(lo,hi,s_update,ng_u,force,ng_o,ave,slx,sly, &
                             slxy,slxx,slyy,uadv,vadv,ng_v,dx,dt,sx,sy,ng_f,sedgex,sedgey,ng_e,bc)

    integer        ,intent(in   ) :: lo(:), hi(:), ng_u, ng_o, ng_f, ng_e, ng_v
    real(kind=dp_t),intent(inout) ::   s_update(lo(1)-ng_u:,lo(2)-ng_u:)
    real(kind=dp_t),intent(in   ) ::      force(lo(1)-ng_o:,lo(2)-ng_o:)
    real(kind=dp_t),intent(in   ) ::  ave(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) ::  slx(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) ::  sly(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) :: slxy(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) :: slxx(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) :: slyy(lo(1)- 1:,lo(2)- 1:)
    real(kind=dp_t),intent(in   ) :: uadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) :: vadv(lo(1)-ng_v:,lo(2)-ng_v:)
    real(kind=dp_t),intent(in   ) :: sx(lo(1)-ng_f:,lo(2)-ng_f:)
    real(kind=dp_t),intent(in   ) :: sy(lo(1)-ng_f:,lo(2)-ng_f:)
    real(kind=dp_t),intent(inout) :: sedgex(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t),intent(inout) :: sedgey(lo(1)-ng_e:,lo(2)-ng_e:)
    real(kind=dp_t),intent(in   ) ::   dx(:),dt
    integer        ,intent(in   ) :: bc(:,:)

    real(kind=dp_t),allocatable ::    gamp(:)
    real(kind=dp_t),allocatable ::    gamm(:)
    real(kind=dp_t),allocatable ::      xm(:)
    real(kind=dp_t),allocatable ::      ym(:)
    real(kind=dp_t),allocatable ::       c(:,:)


    real(kind=dp_t) :: hx,hy,dt3rd,hxs,hys,hxinv,hyinv
    real(kind=dp_t) :: vtrans,stem,vaddif,vdif
    real(kind=dp_t) :: isign, jsign, force_local
    integer i,j,is,ie,js,je
    integer iup,jup

    real(kind=dp_t) :: u1,u2,v1,v2,uu,vv
    real(kind=dp_t) :: eps
    real(kind=dp_t), parameter :: two3rd = 2.d0/3.d0

    allocate( gamp(lo(1)-1:hi(1)+1))
    allocate( gamm(lo(1)-1:hi(1)+1))
    allocate(   xm(lo(1)-1:hi(1)+1))
    allocate(   ym(lo(1)-1:hi(1)+1))
    allocate(    c(lo(1)-1:hi(1)+1,4))

    eps = 1.d-8

    is = lo(1)
    js = lo(2)
    ie = hi(1)
    je = hi(2)
    hx = dx(1)
    hy = dx(2)
    hxinv = 1/hx
    hyinv = 1/hy

    dt3rd = dt / 3.d0

    ! note: the reconstructed function has now the form 
    ! ave + (x-xi)*slx + (y-yj)*sly + (x-xi)(y-yj)*slxy
    ! + (x-xi)^2*slxx + (y-yj)^2*slyy

    do j = js,je 

       ! ******************************* 
       ! calculate Gamma plus for flux F

       do i = is-1,ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j+1)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j
             jsign = 1.d0
             u2 = uadv(i+1,j)
          else 
             jup   = j+1
             jsign = -1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j+1) .gt. 0) then
                u2 = uadv(i+1,j+1)
             end if
          end if

          vv = vadv(iup,j+1)

          hxs = hx*isign
          hys = hy*jsign

          ! quadrature rules:
          ! midpoint for const and linear terms
          ! average of midpoint of edges for bilin and quadr terms

          gamp(i) = ave(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slx(iup,jup) +   &
               (hys*.5 -    vv*dt/3.d0)*sly(iup,jup) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slxy(iup,jup)/12.d0 + &
               ((hxs-dt*u1)**2 + (hxs-dt*u2)**2 + (hxs-dt*(u1+u2))**2)*slxx(iup,jup)/12.d0 + &
               (hys**2+2.*(hys-dt*vv)**2)*slyy(iup,jup)/12.d0

       enddo

       ! end of calculation of Gamma plus for flux F
       ! ****************************************


       ! *****************************************
       ! calculate Gamma minus for flux F

       do i = is-1,ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
          else
             iup   = i+1
             isign = -1.d0
          end if

          vtrans = vadv(iup,j)
          u1 = uadv(i+1,j)
          if (vtrans .gt. 0) then 
             jup   = j-1
             jsign = 1.d0
             u2 = 0.
             if (uadv(i+1,j)*uadv(i+1,j-1) .gt. 0) then
                u2 = uadv(i+1,j-1)
             end if
          else 
             jup   = j
             jsign = -1.d0
             u2 = uadv(i+1,j)
          end if


          vv = vadv(iup,j)

          hxs = hx*isign
          hys = hy*jsign

          gamm(i) = ave(iup,jup)+     &
               (hxs*.5 - (u1+u2)*dt/3.d0)*slx(iup,jup) +   &
               (hys*.5 -    vv*dt/3.d0)*sly(iup,jup) +    &
               (3.*hxs*hys-2.*(u1+u2)*dt*hys-2.*vv*hxs*dt+  &
               vv*(2.*u2+u1)*dt*dt)*slxy(iup,jup)/12.d0 + &
               ((hxs-dt*u1)**2 + (hxs-dt*u2)**2 + (hxs-dt*(u1+u2))**2)*slxx(iup,jup)/12.d0 + &
               (hys**2+2.*(hys-dt*vv)**2)*slyy(iup,jup)/12.d0


       enddo

       ! end of calculation of Gamma minus for flux F
       ! ****************************************


       ! *********************************
       ! calculate sedgex

       do i = is-1, ie 

          if (uadv(i+1,j) .gt. 0) then 
             iup   = i
             isign = 1.d0
             force_local = force(i,j)
          else
             iup   = i+1
             isign = -1.d0
             force_local = force(i+1,j)
          end if

          ! gamm and gamp are updated to handle quadratic terms
          ! stem is updated to handle quadratic terms (using a 1d 2-point Gaussian formula)
          ! this automatically should give `correct' volume integrals (in the non-div free case)

          hxs = hx*isign

          vdif = 0.5d0*dt*(vadv(iup,j+1)*gamp(i) -  &
               vadv(iup,j)*gamm(i) ) * hyinv
          stem = ave(iup,j) + (hxs - uadv(i+1,j)*dt)*0.5d0*slx(iup,j) + &
               0.5d0*slxx(iup,j)*( (hxs/2. - uadv(i+1,j)*dt*(1.+sqrt(3.))/(2.*sqrt(3.)))**2 + &
               (hxs/2. + uadv(i+1,j)*dt*(1.-sqrt(3.))/(2.*sqrt(3.)))**2   ) + &
               slyy(iup,j)*hy*hy/12.d0
          vaddif = stem*0.5d0*dt*(uadv(iup+1,j) - uadv(iup,j)) * hyinv
          sedgex(i+1,j) = stem - vdif - vaddif + 0.5d0*dt*force_local

       enddo
    enddo

    if (bc(1,1) .eq. EXT_DIR) then
       sedgex(lo(1),lo(2):hi(2)) = sx(lo(1),lo(2):hi(2))
    end if

    if (bc(1,2) .eq. EXT_DIR) then
       sedgex(hi(1)+1,lo(2):hi(2)) = sx(hi(1)+1,lo(2):hi(2))
    end if

    ! end of calculation of sedgex
    ! *************************************

    do j = js-1,je 

       ! ********************************** 
       ! calculate Gamma plus for flux G


       do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i+1,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i
             isign = 1.d0
             v2 = vadv(i,j+1)
          else
             iup   = i+1
             isign = -1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i+1,j+1) .gt. 0) then
                v2 = vadv(i+1,j+1)
             end if
          end if


          uu = uadv(i+1,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamp(i) = ave(iup,jup)+ &
               (hys*.5 - (v1+v2)*dt/3.)*sly(iup,jup) +   &
               (hxs*.5 - uu*dt/3.)*slx(iup,jup) + &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slxy(iup,jup)/12.d0 + &
               (hxs**2 + 2.*(hxs-uu*dt)**2)*slxx(iup,jup)/12.d0 + &
               ((hys-v1*dt)**2+(hys-v2*dt)**2+(hys-(v1+v2)*dt)**2)*slyy(iup,jup)/12.d0

       enddo

       ! end of calculation of Gamma plus for flux G
       ! ****************************************


       ! *****************************************
       ! calculate Gamma minus for flux G

       do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
          else
             jup   = j+1
             jsign = -1.d0
          end if

          vtrans = uadv(i,jup)
          v1 = vadv(i,j+1)
          if (vtrans .gt. 0.d0) then
             iup   = i-1
             isign = 1.d0
             v2 = 0.
             if (vadv(i,j+1)*vadv(i-1,j+1) .gt. 0) then
                v2 = vadv(i-1,j+1)
             end if
          else
             iup   = i
             isign = -1.d0
             v2 = vadv(i,j+1)
          end if


          uu = uadv(i,jup)       

          hxs = hx*isign
          hys = hy*jsign

          gamm(i) = ave(iup,jup) +    &
               (hys*.5 - (v1+v2)*dt/3.)*sly(iup,jup) +    &
               (hxs*.5 - uu*dt/3.)*slx(iup,jup) +   &
               (3.*hxs*hys-2.*(v1+v2)*dt*hxs-2.*uu*hys*dt+  &
               (2.*v2+v1)*uu*dt*dt)*slxy(iup,jup)/12.d0 + &
               (hxs**2 + 2.*(hxs-uu*dt)**2)*slxx(iup,jup)/12.d0 + &
               ((hys-v1*dt)**2+(hys-v2*dt)**2+(hys-(v1+v2)*dt)**2)*slyy(iup,jup)/12.d0

       enddo

       ! end of calculation of Gamma minus for flux G
       ! ****************************************


       ! *********************************
       ! calculate sedgey

       do i = is,ie 

          if (vadv(i,j+1) .gt. 0) then 
             jup   = j
             jsign = 1.d0
             force_local = force(i,j)
          else
             jup   = j+1
             jsign = -1.d0
             force_local = force(i,j+1)
          end if

          hys = hy*jsign

          vdif = 0.5d0*dt* &
               (uadv(i+1,jup)*gamp(i)-uadv(i,jup)*gamm(i))*hxinv
          stem = ave(i,jup) + (hys - vadv(i,j+1)*dt)*0.5d0*sly(i,jup) + &
               0.5d0*slyy(i,jup)*( (hys/2. - vadv(i,j+1)*dt*(1.+sqrt(3.))/(2.*sqrt(3.)))**2 + &
               (hys/2. + vadv(i,j+1)*dt*(1.-sqrt(3.))/(2.*sqrt(3.)))**2   ) + &
               slxx(i,jup)*hx*hx/12.d0
          vaddif = stem*0.5d0*dt*(vadv(i,jup+1) - vadv(i,jup))*hyinv
          sedgey(i,j+1) = stem - vdif - vaddif + 0.5d0*dt*force_local

       enddo
    enddo

    if (bc(2,1) .eq. EXT_DIR) then
       sedgey(lo(1):hi(1),lo(2)) = sy(lo(1):hi(1),lo(2))
    end if

    if (bc(2,2) .eq. EXT_DIR) then
       sedgey(lo(1):hi(1),hi(2)+1) = sy(lo(1):hi(1),hi(2)+1)
    end if

    ! end of calculation of sedgey
    ! *************************************

    deallocate(gamp,gamm,xm,ym,c)

  end subroutine bdsconc_quad_2d

end module bds_module
