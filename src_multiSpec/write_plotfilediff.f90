module write_plotfilediff_module

  use ml_layout_module
  use multifab_module
  use fabio_module
  use probin_multispecies_module, only: nspecies
  use probin_common_module, only: prob_lo, prob_hi

  implicit none

contains
  
  subroutine write_plotfilediff(mla,name,rho,istep,dx,time)

    type(ml_layout),    intent(in)  :: mla
    character(len=*),   intent(in)  :: name
    type(multifab),     intent(in)  :: rho(:)
    integer,            intent(in)  :: istep
    real(kind=dp_t),    intent(in)  :: dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plotfile_name
    integer                         :: n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"write_plotfilediff")

    nlevs = mla%nlevel
  
    allocate(plot_names(nspecies))
    allocate(plotdata(nlevs))
 
    ! write density with species index
    do n=1,nspecies
       write(plot_names(n),'(a,i0)') "rho", n
    enddo

    ! build plotdata for nspecies and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),nspecies,0)
    enddo
    
    ! copy the state into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1,rho(n),1,nspecies,0)
    enddo
    
    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='(a,i8.8)') name, istep
    !if ( parallel_IOProcessor() ) then
    !  write(*,'(2A)') "Saving PLOT FILEs to directory ", trim(plotfile_name)
    !  write(*,*)
    !end if
    
    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), trim(plotfile_name), & 
                                   plot_names, mla%mba%pd(1), prob_lo, prob_hi, & 
                                   time, dx(1,:))

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    enddo
    
    deallocate(plotdata)
    deallocate(plot_names)

    call destroy(bpt)

  end subroutine write_plotfilediff
  
  subroutine write_plotfile_onecomp(mla,name,Temp,istep,dx,time)

    type(ml_layout),    intent(in)  :: mla
    character(len=*),   intent(in)  :: name
    type(multifab),     intent(in)  :: Temp(:)
    integer,            intent(in)  :: istep
    real(kind=dp_t),    intent(in)  :: dx(:,:),time

    ! local variables
    character(len=20), allocatable  :: plot_names(:)
    character(len=20)               :: plotfile_name
    integer                         :: n,nlevs

    ! multifab of size nlevs  
    type(multifab), allocatable     :: plotdata(:)

    type(bl_prof_timer), save :: bpt

    call build(bpt,"write_plotfile_onecomp")

    nlevs = mla%nlevel
  
    allocate(plot_names(1))
    allocate(plotdata(nlevs))
 
    ! write density with species index
    write(plot_names(1),'(a,i0)') "Scalar"

    ! build plotdata for nspecies and 0 ghost cells
    do n=1,nlevs
       call multifab_build(plotdata(n),mla%la(n),1,0)
    enddo
    
    ! copy the state into plotdata
    do n = 1,nlevs
       call multifab_copy_c(plotdata(n),1,Temp(n),1,1,0)
    enddo
    
    ! define the name of the plotfile that will be written
    write(unit=plotfile_name,fmt='(a,i8.8)') name, istep
        
    ! write the plotfile
    call fabio_ml_multifab_write_d(plotdata, mla%mba%rr(:,1), trim(plotfile_name), & 
                                   plot_names, mla%mba%pd(1), prob_lo, prob_hi, & 
                                   time, dx(1,:))

    ! make sure to destroy the multifab or you'll leak memory
    do n=1,nlevs
       call multifab_destroy(plotdata(n))
    enddo
    
    deallocate(plotdata)
    deallocate(plot_names)

    call destroy(bpt)

  end subroutine write_plotfile_onecomp

end module write_plotfilediff_module
