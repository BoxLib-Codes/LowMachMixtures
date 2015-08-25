module ParallelRNGs
   ! Here BoxLib is used to seed the RNGs
   ! Building this code requires having BoxLib installed
   
   use BoxLibRNGs
   use Random_Numbers ! Used to generate unpredictable seeds or initial seeds for Marsenne twister
   use parallel
   use bl_error_module
   use vector_i_module
   
   implicit none
   private
   
   public :: SeedParallelRNG
   
contains 

  ! Note: boxlib_initialize must  have already been called!
  subroutine SeedParallelRNG (rootseed)
    integer, intent(inout) :: rootseed

    integer :: iseed, i, ourseed(1), NProcs
    integer, allocatable :: seeds(:)
    type(vector_i) vseeds

    NProcs = parallel_nprocs()
    allocate(seeds(NProcs))
    
    call build(vseeds)

    if ( parallel_ioprocessor() ) then
       !
       ! We'll get NProcs unique integers with which to seed MT.
       !
       if ( rootseed == 0 ) then
          call UnpredictableSeeds(rootseed)
       else
          call RandomSeeds(rootseed)
       end if
       write(*,*) "Using initial seed on proc(0): ROOT_SEED=", rootseed

       do
          call RandomUniform(iseed)

          if ( iseed == 0 ) cycle

          do i = 1, size(vseeds)
             if ( at(vseeds,i) == iseed ) exit
          end do

          if ( i > size(vseeds) ) then
             call push_back(vseeds,iseed)
          end if

          if ( size(vseeds) == NProcs ) exit
       end do

       call bl_assert(size(vseeds) == NProcs, 'should have NProcs seeds')

       do i = 1, NProcs
          seeds(i) = at(vseeds,i)
       end do
    end if

    call parallel_scatter(seeds, ourseed, 1)
    !
    ! Now set'm ...
    !
    call SeedRNG(ourseed(1))

    if ( parallel_ioprocessor() ) then

       print*, ' '
       if ( NProcs > 1 ) then
          do i = 1,NProcs
             print*, "CPU# ", (i-1), ", SEED = ", seeds(i)
          end do
       else
          print*, "Initial SEED = ", ourseed(1)
       end if
       print*, ' '

    end if

    call destroy(vseeds)

  end subroutine
  
end module
