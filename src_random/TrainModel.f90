program TrainModel
   use Precision
   use BoxLibRNGs
   use HydroGridModule
   implicit none
   
   integer, parameter :: wp=r_dp, ip=i_sp ! Working precision

   integer(i_dp) :: step, n_steps, n_steps_snapshots, n_steps_save, n_skip_steps
   
   integer :: seed
   integer(ip) :: n_tracks, track, direction, n_passrs_per_train, n_passengers, max_density
   ! Mass of each passenger is assumed m=1, and dx=1, so density is an integer
   integer(ip), dimension(:), allocatable :: density
   real(wp), dimension(:), allocatable :: momentum, massflux, momentumflux
   real(wp) :: v_platform, time, dice, dp, last_sample_time, dn_platform
   
   type(HydroGrid) :: grid
   integer :: nameListFile=173
   
   !----------------------------------

   read(*,*) n_tracks, n_steps, seed
   read(*,*) n_passrs_per_train, v_platform, dn_platform
   read(*,*) n_steps_snapshots, n_steps_save, n_skip_steps
   
   call SeedRNG(seed)
   
   allocate(density(0:n_tracks+1), momentum(0:n_tracks+1))
   allocate(massflux(n_tracks+1), momentumflux(n_tracks+1))

   ! The average time interval per step is 1/n_passengers
   !open (nameListFile, file = "TrainModel.nml", status="old", action="read")
   nameListFile=5
   call createHydroAnalysis (grid, nCells=(/1,n_tracks,1/), nSpecies=1, nVelocityDimensions=2, &
      isSingleFluid=.true., systemLength=(/1.0_wp,1.0_wp*n_tracks,1.0_wp/), &
      heatCapacity=(/1.0_wp/), timestep=n_steps_snapshots*(1.0_wp/n_passrs_per_train)/n_tracks, &
      fileUnit=nameListFile, structFactMultiplier=real(n_passrs_per_train,wp))
   !close(nameListFile)      
   
   ! Initialize with shear velocity profile
   do track=0, n_tracks+1
      density(track)=nint(n_passrs_per_train+track*dn_platform/(n_tracks+1))
      momentum(track)=density(track)*track*v_platform/(n_tracks+1)
   end do
   ! We will keep track of the total momentum imparted to the platforms
   momentum(0)=0
   momentum(n_tracks+1)=0
   
   n_passengers=sum(density) ! Total number of passengers
   max_density=maxval(density) ! For rejection MC
   massflux=0 ! Keep track of f, the mass flux of passengers between trains
   momentumflux=0 ! Keep track of "v*f", the momentum flux between trains
   time=0
   last_sample_time=0
   step=0
   do
      ! No need for Poisson increments for steady state
      time = time + 1.0_wp/max_density/(n_tracks+1)
         ! Because we use rejection here, the number of passengers used
         ! is the ficticious one, assuming density=max_density
         ! The end platforms are selected at half the rate, so add +1
      
      ! Choose the track of the next passenger to jump
      call UniformRNG(dice)
      ! This chooses the platforms at half the rate:
      track = nint(dice*(n_tracks+1)) ! Nearest integer
      
      ! Accept or reject this track
      call UniformRNG(dice)
      if(density(track)<=dice*max_density) cycle ! Reject      
            
      ! dp is the momentum of the jumping passenger
      ! mass of passanger is unity here      
      if(track==0) then
         n_passengers=n_passengers+1

         density(1) = density(1) + 1
         max_density=max(max_density, density(1))

         dp = 0 ! Platform has zero velocity
         momentum(0) = momentum(0) - dp
         momentum(1) = momentum(1) + dp
         
         massflux(1) = massflux(1) + 1
         momentumflux(1) = momentumflux(1) + dp

      else if(track==n_tracks+1) then
         n_passengers=n_passengers+1

         density(n_tracks) = density(n_tracks) + 1
         max_density=max(max_density, density(n_tracks))

         dp = v_platform ! Plaform has fixed velocity
         momentum(n_tracks+1) = momentum(n_tracks+1) - dp
         momentum(n_tracks) = momentum(n_tracks) + dp
         
         massflux(n_tracks+1) = massflux(n_tracks+1)  - 1
         momentumflux(n_tracks+1) = momentumflux(n_tracks+1)  - dp

      else

         ! No danger of division by zero here due to rejection MC
         dp = momentum(track)/density(track) ! Momentum of jumping passenger

         ! Choose a direction to jump in
         call UniformRNG(dice)
         if(dice<0.5_wp) then
            direction=-1 ! Jump to the left
            massflux(track) = massflux(track) - 1
            momentumflux(track) = momentumflux(track) - dp
         else
            direction=1 ! Jump to the right
            massflux(track+1) = massflux(track+1) + 1
            momentumflux(track+1) = momentumflux(track+1) + dp
         end if      
         
         density(track) = density(track) - 1
         momentum(track) = momentum(track) - dp
         momentum(track+direction) = momentum(track+direction) + dp                  
               
         if((track+direction)==0) then
            n_passengers=n_passengers-1
         else if((track+direction)==n_tracks+1) then
            n_passengers=n_passengers-1
         else
            density(track+direction) = density(track+direction) + 1
            max_density=max(max_density, density(track+direction))            
         end if
         
      end if      
      !if(n_passengers /= sum(density)) stop "Imbalanced n_passengers"
      
      step=step+1 ! We count actual hops here
      if(step<n_skip_steps) then
         cycle
      else if(step==n_skip_steps) then
         step=0
         n_skip_steps=-1
      end if
            
      if(mod(step,n_steps_snapshots)==0) then
         ! Verify conservation
         if(mod(step,100000*n_steps_snapshots)==0) then
            !write(*,*) step, n_passengers, sum(momentum)
            write(*,*) step/n_steps_snapshots, (/momentum(0),-momentum(n_tracks+1)/) / &
               (time*v_platform/(n_tracks+1))/n_passrs_per_train
         end if      
         
         ! Calculate average fluxes (transport per unit time)
         ! But these are assumed to be white in time, so divide by sqrt(dt)
         massflux=massflux/sqrt(time-last_sample_time)
         ! We treat this one as non-white
         momentumflux=momentumflux/(time-last_sample_time)

         ! Do some analysis using the HydroGrid code
         ! We pretend here that the massflux f is a face-centered (staggered) velocity (vy)
         call updateHydroAnalysisStaggered (grid, nGhost=0, nGhostScalar=0, &
                  density=real(density(1:n_tracks),wp), &
                  vx=momentum(1:n_tracks)/max(1,density(1:n_tracks)), vy=massflux(1:n_tracks))
         
         ! Reset some variables:
         max_density=maxval(density)
         last_sample_time=time
         massflux=0
         momentumflux=0
         
      end if
      
      if((n_steps_save>0).and.mod(step,n_steps_save)==0) then
         call writeToFiles(grid, id=int(step/n_steps_save,ip))
      end if   
      
      if(step>=n_steps) exit
   end do
   
   write(*,*) "Averaged over ", step/n_steps_snapshots, " samples"
   call writeToFiles(grid)
   call destroyHydroAnalysis(grid)

   if(.false.) then
      do track=0, n_tracks+1
         write(11,*) density(track)
      end do
      do track=1, n_tracks+1
         write(12,*) massflux(track)/time
      end do   
      do track=1, n_tracks+1
         write(13,*) momentumflux(track)/time
      end do
   end if
      
   ! The effective shear rate is v_platform/(n_tracks+1)
   write(0,*) (/momentum(0),-momentum(n_tracks+1)/) / &
      (time*v_platform/(n_tracks+1))/n_passrs_per_train      
         
end program
