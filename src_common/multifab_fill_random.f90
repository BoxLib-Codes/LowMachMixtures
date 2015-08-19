module multifab_fill_random_module

  use multifab_module
  use BoxLibRNGs
  use bl_error_module

  implicit none

  private

  public :: multifab_fill_random

contains

  ! fill a multifab with random numbers
  subroutine multifab_fill_random(mfab, comp, variance, variance_mfab, variance_mfab2)
    type(multifab), intent(inout)           :: mfab(:)
    integer       , intent(in   ), optional :: comp ! Only one component
    real(dp_t)    , intent(in   ), optional :: variance
    type(multifab), intent(in   ), optional :: variance_mfab(:)
    type(multifab), intent(in   ), optional :: variance_mfab2(:)

    integer :: n,box

    real(kind=dp_t), pointer :: fp(:,:,:,:), fpvar(:,:,:,:)

    !--------------------------------------
    do n=1,size(mfab)
       do box = 1, nfabs(mfab(n))
          if(present(comp)) then
             fp => dataptr(mfab(n),box,comp,1)
          else
             fp => dataptr(mfab(n),box)
          end if

          call NormalRNGs(fp, size(fp)) ! Fill the whole grid with random numbers

          if(present(variance_mfab)) then 
             ! Must have same distribution
             if (mfab(n)%ng .ne. variance_mfab(n)%ng) then
                call bl_error('ERROR: multifab_fill_random; mfab and variance_mfab require same # of ghost cells')
             end if
             fpvar => dataptr(variance_mfab(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if

          if(present(variance_mfab2)) then 
             ! Must have same distribution
             if (mfab(n)%ng .ne. variance_mfab2(n)%ng) then
                call bl_error('ERROR: multifab_fill_random; mfab and variance_mfab2 require same # of ghost cells')
             end if
             fpvar => dataptr(variance_mfab2(n),box,1,size(fp,4))
             fp=sqrt(fpvar)*fp
          end if

          if(present(variance)) then
             fp=sqrt(variance)*fp
          end if

       end do
    end do

  end subroutine multifab_fill_random

end module multifab_fill_random_module
