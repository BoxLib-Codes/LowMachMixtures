module multifab_filter_module

  use multifab_module

  implicit none

  private

  public :: multifab_filter

contains

  subroutine multifab_filter(mfab, width)

    ! Note: This routine assumes ghost values are valid on entry but does NOT 
    !       sync ghost values at the end.
    !       If ghost values need to remain consistent do a multifab_fill_boundary 
    !       after calling multifab_filter!

    type(multifab), intent(inout) :: mfab
    integer       , intent(in   ) :: width

    ! local
    integer :: i,j,k,box,dm
    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    integer :: lo(3), hi(3), lo_g(3), hi_g(3)
    
    type(bl_prof_timer),save :: bpt

    call build(bpt,"multifab_filter")

    lo=1
    hi=1
    lo_g=1
    hi_g=1
  
    dm = multifab_get_dim(mfab)

    !--------------------------------------
    FilterX: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       !write(*,*) "lo=", lo, " hi=", hi
       
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       
       !write(*,*) "lo_g=", lo_g, " hi_g=", hi_g

       if(.not.(all(lo_g(1:dm)<=lo(1:dm)-width) .and. all(hi_g(1:dm)>=hi(1:dm)+width) ) ) then
         call bl_error("Filtering requires at least width cells!")
       end if
       
       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array
              
       ! Filter along x:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i-1,j,k,:)/4+fp(i+1,j,k,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                     -fp(i-2,j,k,:)/16+fp(i-1,j,k,:)/4 &
                     -fp(i+2,j,k,:)/16+fp(i+1,j,k,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i-1,j,k,:)/32 - 7*fp(i-2,j,k,:)/64 + fp(i-3,j,k,:)/32 - fp(i-4,j,k,:)/256 &
                     + 7*fp(i+1,j,k,:)/32 - 7*fp(i+2,j,k,:)/64 + fp(i+3,j,k,:)/32 - fp(i+4,j,k,:)/256 )
                case default
                  call bl_error("Input width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterX
   
    !--------------------------------------
    ! Update ghost values to reflect filtering
    call multifab_fill_boundary(mfab)
   
    FilterY: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       

       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array

       ! Filter along y:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i,j-1,k,:)/4+fp(i,j+1,k,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                    -fp(i,j-2,k,:)/16+fp(i,j-1,k,:)/4 &
                    -fp(i,j+2,k,:)/16+fp(i,j+1,k,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i,j-1,k,:)/32 - 7*fp(i,j-2,k,:)/64 + fp(i,j-3,k,:)/32 - fp(i,j-4,k,:)/256 &
                     + 7*fp(i,j+1,k,:)/32 - 7*fp(i,j+2,k,:)/64 + fp(i,j+3,k,:)/32 - fp(i,j+4,k,:)/256 )
                case default
                  call bl_error("Input width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterY
    !--------------------------------------

    if(dm<=2) return ! We are done
    
    !--------------------------------------
    call multifab_fill_boundary(mfab) ! Update ghost values to reflect filtering
   
    FilterZ: do box = 1, nfabs(mfab)

       fp => dataptr(mfab,box) ! Including ghosts
       
       ! Without ghosts:
       lo(1:dm) = lwb(get_ibox(mfab,box))
       hi(1:dm) = upb(get_ibox(mfab,box))
       ! With ghosts:
       lo_g(1:dm) = lwb(get_pbox(mfab,box))
       hi_g(1:dm) = upb(get_pbox(mfab,box))       

       allocate(fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3), &
                      lbound(fp,4):ubound(fp,4))) ! Temporary array

       ! Filter along z:
       do k = lo(3),hi(3)
          do j = lo(2),hi(2)
             do i = lo(1),hi(1)
                select case(width)
                case(1)
                   fpvar(i,j,k,:) = ( fp(i,j,k,:)/2 &
                     +fp(i,j,k-1,:)/4+fp(i,j,k+1,:)/4 )
                case(2)
                   fpvar(i,j,k,:) = ( 5*fp(i,j,k,:)/8 &
                    -fp(i,j,k-2,:)/16+fp(i,j,k-1,:)/4 &
                    -fp(i,j,k+2,:)/16+fp(i,j,k+1,:)/4 )
                case(4)
                   fpvar(i,j,k,:) = ( 93*fp(i,j,k,:)/128 &
                     + 7*fp(i,j,k-1,:)/32 - 7*fp(i,j,k-2,:)/64 + fp(i,j,k-3,:)/32 - fp(i,j,k-4,:)/256 &
                     + 7*fp(i,j,k+1,:)/32 - 7*fp(i,j,k+2,:)/64 + fp(i,j,k+3,:)/32 - fp(i,j,k+4,:)/256 )
                case default
                  call bl_error("Input width must be 0,1,2 or 4")
                end select      
             end do
          end do
       end do
       fp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = fpvar(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:)

       deallocate(fpvar)

    end do FilterZ
    !--------------------------------------

    call destroy(bpt)

  end subroutine multifab_filter

end module multifab_filter_module
