program ViscousBurgers
   ! This solves the viscous Burgers or advection-diffusion equation:
   ! u_t + v*u_x = D*u_xx + [epsilon*2*D*A(u)*W]_x
   ! where v is an advection speed:
   ! v=c*u if c>0 (Burgers), or v=abs(c) if c<=0 (advection-diffusion)

   ! In the code c=advspeed, D=diffcoeff, epsilon=varsc
   
   ! If this is EWR Burgers, then u is in [-1,1] and the actual density is
   ! rho = (1-u)/2 and is in [0,1]
   
   ! And we have the following cases:
   ! noise = (note that physical prefactors such as 2 are omitted here)
   ! If negative, the noise is calculated at the initial condition only and kept fixed
   !       0 => A(u)=1 ! Additive noise
   !       1 => A(u)=u ! As in simple diffusion
   !       2 => A(u)=(1-u)/2 ! Corresponding to A=rho for EWR
   !       3 => A(u)=u*(1-u) ! As in concentration equation
   !       4 => A(u)=(1-u^2)/4 ! ERW, corresponding to A=rho*(1-rho)
   !       5 => A(u)=u^2 ! As in temperature equation

   use Precision
   use BoxLibRNGs
   use HydroGridModule
   use TridiagonalSolvers
   implicit none
   
   integer, parameter :: wp=r_dp, ip=i_sp ! Working precision

   integer(i_dp) :: step, n_steps, n_steps_snapshots, n_steps_save, n_skip_steps
   
   integer :: seed=0, n_runs, run
   integer :: n_cells=100, cell, face
   real(wp), dimension(:), allocatable :: density, ddensity, old_density, &
      correlation, mean, massflux, stochforcing, denscopy, ddensity_copy
   real(wp), dimension(:,:), allocatable :: rngs
   real(wp) :: time, dt=1.0, dx=1.0, fdensity, diffcoeff=1.0, advspeed=0
   real(wp) :: bcs(2)=0.5_wp ! Values at boundaries
   real(wp) :: varsc=1, gamma_sign=1, force_const=1
   integer :: noise=0, initcond=0 ! Initial condition:
      ! 0=constant, 1=linear, 2=moving shock, 3=stationary shock, 4=periodic solution
      ! Negative adds a suspended particle: -1=free particle, -2=const. forcing, -3=Harmonically bound
   integer :: kernel_width=3 ! Can be a multiple of 3 or 4    
   integer :: integrator=1 ! 0=Euler
      ! Explicit: 1=Trapezoidal PC, 2=Midpoint PC, 3=RK3, 4=old RK3
      ! Semi-implicit: 10=Euler, 11=CN trapezoidal, 12=CN midpoint,
      ! 13=L-stable trapezoidal, 14=L-stable midpoint, 15=single CN step
      ! If negative a single random number per step is used (same W in all stages)
   integer :: particle_integrator=-1
      ! 0=Lie splitting, -1=Strang trapezoidal, -2=Strang midpoint, -3=Strang Euler,
      ! With fluid: 1=Explicit, 2=Implicit (not implemented!),
      ! 10=Staggered Euler, 11=Staggered trap, 12=Staggered mid
      ! Boyce's explicit predictor-corrector: 21=Trapezoidal, 22=Midpoint
   logical :: periodic=.false.
   integer :: interpolations(2)=0 ! How to interpolate from cells to faces
      ! interpolations(1) for advection speed, =0 for average, =1 for Hamiltonian
      ! interpolations(2) for noise,
      !    =0 for average of state, =1 for average of variances, >=2 for fancy stuff
   
   integer :: pint_copy, int_copy, n_init_steps
   real(wp) :: l2r, r2l, adv_vel, r_tmp, step_length, pi, system_length, dt_copy
   real(wp) :: mincenter, shockpos ! Positions of shock
   real(wp) :: w1, w2, w3, w4, w5, w6 ! Weights
   
   ! Particle data
   real(wp) :: initial_position, position, position_next, periodic_shift, poscopy, &
      velocity, velocity_next, local_v, force, pos_stag_last, pos_stag_next, &
      dposition, dposcopy, forcecopy, position_shift, shifted_position
   
   type(HydroGrid) :: grid
   integer :: nameListFile=173
   
   namelist /Burgers/ n_cells, periodic, seed, & ! System size
      advspeed, diffcoeff, noise, varsc, & ! Equation form
      initcond, bcs, & ! Values used for BCs
      dx, dt, interpolations, integrator, & ! Discretization
      n_steps, n_steps_snapshots, n_steps_save, n_skip_steps, & ! Iteration counts
      n_runs, gamma_sign, kernel_width, particle_integrator, force_const

   pi=4*atan(1.0_wp)
      
   
   ! ===================== UNFINISHED =========================
   !open (nameListFile, file = "ViscousBurgers.nml", status="old", action="read")
   nameListFile=5 ! Standard input (stdin) 
   
   read(unit=nameListFile,nml=Burgers)
   
   call SeedRNG(seed)
   
   allocate(density(0:n_cells+1), denscopy(0:n_cells+1), old_density(0:n_cells+1))
   allocate(ddensity(1:n_cells), ddensity_copy(1:n_cells))
   allocate(mean(0:n_cells+1), correlation(0:n_cells+1))
   allocate(massflux(1:n_cells+1), stochforcing(1:n_cells+1))
   allocate(rngs(1:n_cells+1,2))

   if(n_steps_snapshots>0) then
      call createHydroAnalysis (grid, nCells=(/1,n_cells,1/), nSpecies=1, nVelocityDimensions=2, &
         isSingleFluid=.true., systemLength=(/1.0_r_wp,dx*n_cells,1.0_r_wp/), &
         heatCapacity=(/1.0_wp/), timestep=n_steps_snapshots*dt, &
         fileUnit=nameListFile, structFactMultiplier=1.0/varsc)
   end if   
      
   !close(nameListFile)      
   
   if(dx<0) then
      if(initcond==4) then
         dx=2*pi/n_cells
      else
         dx=abs(dx)/n_cells
      end if   
   end if
   
   if(dt<0) then ! For convergence studies
      dt=abs(dt)/n_cells
      write(*,*) "Using timestep dt=", dt
   end if

   ! Set various weights required by the algorithms:
   select case(abs(integrator))
   case(3) ! RK3 new weights better
      w1=(2*sqrt(2.0_wp)+sqrt(3.0_wp))/5
      w2=(-4*sqrt(2.0_wp)+3*sqrt(3.0_wp))/5
      w3=(sqrt(2.0_wp)-2*sqrt(3.0_wp))/10
      w4=0; w5=0; w6=0;
   case(5) ! RK3 new weights but worse
      w1=(2*sqrt(2.0_wp)-sqrt(3.0_wp))/5
      w2=(-4*sqrt(2.0_wp)-3*sqrt(3.0_wp))/5
      w3=(sqrt(2.0_wp)+2*sqrt(3.0_wp))/10
      w4=0; w5=0; w6=0;
   case(4) ! Old weights (first order only)
      w1=-sqrt(3.0_wp)
      w2=sqrt(3.0_wp)
      w3=0
      w4=0; w5=0; w6=0;
   case(13,14) ! L0-stable
      if(integrator==13) then ! Trapezoidal
         w2=1.0_wp
      else
         w2=0.5_wp
      end if
      ! Possible values of interest are sign=-1,1
      w4=1+gamma_sign*sqrt(2.0_wp)/2
      w1=(0.5_wp-w4)/(1-w4)
      w3=(0.5_wp-w4)/w2
      w5=0.5_wp/w2
      w6=0
   end select
   write(*,*) "Using weights ", w1, w2, w3, w4, w5, w6

RunLoop: do run=1, n_runs ! We can do many simulations in sequence
   ! Set parameters
   if(advspeed<0) adv_vel=abs(advspeed) ! Constant advection velocity

   time=0 ! Initial time
   call InitializeFluid()

   if(run==1) then
      if(advspeed<=0) then
         write(*,*) "ADVECTIVE CFL=", dt*abs(advspeed)/dx, " REYNOLDS=", (abs(advspeed)/dx)/(diffcoeff/dx**2)
      else
         r_tmp=advspeed*maxval(abs(density))/dx
         write(*,*) "ADVECTIVE CFL=", dt*r_tmp, " REYNOLDS=", r_tmp/(diffcoeff/dx**2)
      end if
      write(*,*) "DIFFUSIVE CFL=", dt*diffcoeff/dx**2
      if(initcond>=2) then
         write(*,*) "Shock to move half domain in ", &
            nint( ((n_cells+1)*dx/abs(advspeed)/abs(bcs(1)+bcs(2)))/dt ), " steps"
      end if
   end if
      
   if(initcond<0) then ! Advect a particle around, starting in the middle
      if(periodic) then
         position=system_length/2
      else   
         position = system_length/9 ! system_length/17 ! 2.0_wp ! 2*dx ! 0.5_wp*n_cells*dx
      end if   
      initial_position = position
      periodic_shift = 0 ! For counting periodic boundary crossings
      
      ! For staggered schemes we need to obtain consistent initial conditions
      if(particle_integrator>=10) then ! Stagger positions
      if(.true.) then ! Do simple initialization
      
         call move_particle(dt/2, preliminary=.true.) ! Simple Euler step
         
         ! We need q(1/2) to start the algorithm:
         pos_stag_next = position ! q(1/2)
         position = initial_position ! q(0)
      
      else ! Do more correct initialization by running a more accurate integrator         
      
         int_copy=integrator
         pint_copy=particle_integrator
         dt_copy=dt
         
         ! Use a third-order accurate explicit method to go back by -dt/2:
         dt=min(dt_copy/20, 0.1_wp*dx**2/diffcoeff)
         n_init_steps=nint(dt_copy/2/dt)
         write(*,*) "Running BACKWARD in time with ", n_init_steps, " steps of RK3"

         dt=dt_copy/2/n_init_steps ! Go back in time
         integrator=3
         particle_integrator=1
         
         call InitializeTimestepping()
         do step=1, n_init_steps
            call TakeTimestep()
            !write(31,*) real(time), real(position), real(velocity)
         end do         
         pos_stag_next = position ! q(1/2)

         time=0
         call InitializeFluid()
         position=initial_position
         periodic_shift = 0
         integrator=int_copy
         particle_integrator=pint_copy
         dt=dt_copy 
      
      end if
      end if
      
   end if   
      
   !=====================
   ! Begin the run

   call InitializeTimestepping()
   
   step=0
   correlation=0
   mean=0
   position_shift=periodic_shift
   TimeStepLoop: do
      call TakeTimestep()
      
      step=step+1 ! We count actual hops here
      if(step<n_skip_steps) then
         cycle
      else if(step==n_skip_steps) then
         step=1
         n_skip_steps=-1
      end if

      if(initcond<0) then ! Analyze particle data
         ! The definition of local velocity is unique with no inertia:
         call interpolate(velocity, density, position)
         
         ! Convergence tests of pure advection:
         if(.false.) then ! Linear velocity profile
            
            r_tmp = (bcs(2)-bcs(1))/((n_cells+1)*dx) ! Slope s
            velocity = bcs(1) + initial_position*r_tmp
            write(31,*) real(time), real(position), abs(position - (initial_position &
               + velocity*(exp(time*r_tmp)-1)/r_tmp))
               
         else if(.false.) then ! Quadratic velocity profile
         
            r_tmp = (bcs(2)-bcs(1))/((n_cells+1)*dx)**2 ! Slope s
            velocity = bcs(1) + (bcs(2)-bcs(1))*(initial_position/((n_cells+1)*dx))**2
            write(31,*) real(time), real(position), abs(position - &
               tan(time*sqrt(r_tmp*velocity)+ &
               atan(r_tmp*initial_position/sqrt(r_tmp*velocity)))*sqrt(r_tmp*velocity)/r_tmp)
               
         else ! Non-analytical convergence tests:
            write(0,*) real(time), real(shifted_position), real(velocity)
         end if    
            
      end if

      if(.false.) then ! This is now handled by the HydroGrid code
         correlation=correlation+density*density ! Variance
         !correlation=correlation+density*density(n_cells/2) ! Correlation
         mean=mean+density
      end if   

      ! Verify conservation and other analysis
      if(mod(step,n_steps/8)==0) then
         if(initcond==3) then ! Stationary shock
         
            if(run==1) then
               write(100+step/(n_steps/10),'(A,g17.9)') "# t=", time
            end if   
            !call write_shock()
            call fit_shock(mincenter, shockpos)
            !write(*,*) time," Estimated shock positions: ", real(shockpos), real(mincenter)
            write(100+step/(n_steps/10),'(100(g17.9))') mincenter, shockpos
            
         else if(initcond==2) then ! Moving shock
         
            call write_shock()
            call fit_shock(mincenter, shockpos)
            write(*,*) time," Estimated shock positions: ", real(shockpos), real(mincenter)
            
         else if(initcond==4) then ! Exact periodic solution 
         
            !call periodic_solution_Fay(denscopy,time) ! Triangle wave
            call periodic_solution(denscopy,time) ! A simpler one-mode form
            call write_solution(denscopy)
            ! Write the L1, L2, and Linf norms to a file
            write(21,'(100(g17.9))') time, &
               real(sum(abs(density(1:n_cells)-denscopy(1:n_cells)))/n_cells), &
               real(sqrt(sum((density(1:n_cells)-denscopy(1:n_cells))**2)/n_cells)), &
               real(maxval(abs(density(1:n_cells)-denscopy(1:n_cells))))
            
            write(*,*) step, real(time), " <rho>=", sum(density(1:n_cells))/n_cells
            
         else if(.false.) then ! Move the particle artificially to sample J*L^(-1)*S
            call write_solution(density)
            
            call interpolate(velocity, density, position)
            write(103,*) position/system_length, velocity
            write(*,*) step, real(time), " <rho>=", sum(density(1:n_cells))/n_cells
            if(periodic) then
               position = position + dx/9
            else   
               position = position + system_length/9
            end if
            call fix_particle()

         else
             
            !call write_solution(density)
            write(*,*) step, real(time), " <rho>=", sum(density(1:n_cells))/n_cells
         end if
         ! This is verbose, so use only for debugging stability limits etc.
         !write(*,*) step, real(time), " <rho>=", sum(density(1:n_cells))/n_cells
      end if      
      
      if((n_steps_snapshots>0).and.(mod(step,n_steps_snapshots)==0)) then         
         ! Do some analysis using the HydroGrid code
         ! We pretend here that the massflux f is a face-centered (staggered) velocity (vy)
         call updateHydroAnalysisStaggered (grid, nGhost=0, nGhostScalar=0, &
                  density=density(1:n_cells), &
                  vx=density(1:n_cells), vy=massflux(1:n_cells))                  
      end if
      
      if((n_steps_save>0).and.mod(step,n_steps_save)==0) then
         call writeToFiles(grid, id=int(step/n_steps_save,ip))
      end if   
      
      if(step>=n_steps) exit
   end do TimeStepLoop
   
   if(initcond==2) call write_shock()
   !call write_fitted_shock(mincenter, shockpos)
   !write(*,*) time," Estimated shock positions: ", real(shockpos), real(mincenter)
   
   if(.false.) then
      mean = mean/step   
      correlation = correlation/step - mean*mean ! Variance
      !correlation = correlation/step - mean*mean(n_cells/2) ! Correlation
      do cell=0,n_cells+1
         if(cell/=(n_cells/2)) write(11,*) cell, correlation(cell)
         write(12,*) cell, mean(cell)
         write(13,*) cell, density(cell)
      end do
   end if   
   
   if (n_steps_snapshots>0) then
      if(n_runs>1) then
         call writeToFiles(grid, id=run)
      else   
         write(*,*) "Averaged over ", step/n_steps_snapshots, " samples"
         call writeToFiles(grid)
      end if   
   end if   

   if((run==1).or.(mod(run,max(1,n_runs/10))==0)) then
      write(*,*) "Finished run ", run
   end if   
end do RunLoop

   if (n_steps_snapshots>0) call destroyHydroAnalysis(grid)

contains

subroutine InitializeTimestepping()
   if(noise<=0) then
      ! Initialize the stochastic forcing amplitude
      do face=1, n_cells+1
         select case(interpolations(2))
         case(0)
           fdensity = 0.5_wp*(density(face-1)+density(face)) ! Face-interpolated density
           stochforcing(face) = sqrt(2*A_u(fdensity))
         case(1)  
           stochforcing(face) = sqrt((A_u(density(face-1))+A_u(density(face))))
         case(2)
            if(abs(noise)/=3) stop "Fancy interpolation only implemented for noise=3"
            stochforcing(face) = sqrt((density(face-1)*(1-density(face)) + &
                                       density(face)*(1-density(face-1))))
         end select
      end do
   end if   
end subroutine

subroutine TakeTimestep() ! Main routine in the algorithm

   time = time + dt

   ! Strang-splitting algorithm
   if(initcond<0) then
      poscopy=position ! Save the old position

      select case(particle_integrator)
      case(0) ! Use Lie splitting
         ! Move the particle at the end of the step
      case(-1,-2,-3) ! Use Strang splitting
         call move_particle(dt/2)
      case(21) ! Boyce Trapezoidal   
         ! We need to evaluate the force first
      case(22) ! Boyce Midpoint: Take a predictor step to the middle       
         call move_particle(dt/2, preliminary=.true.)
      case(10,11,12) ! Staggered Euler/Midpoint
         ! Here: pos_stag_last=q(n-1/2), pos_stag_next=q(n+1/2)
         ! We keep the position at midpoint q(n+1/2) during the fluid update:
         position = pos_stag_next ! q(n+1/2)
      case(1) ! No splitting
         ! Position will be updated together with velocity
      case(2) ! Implici -- not yet implemented
         ! Do not move the particle at all
      case default
         stop "Unimplemented choice for particle_integrator"   
      end select

      call calc_force(position, force)

      if(particle_integrator==21) then
         forcecopy=force
         call move_particle(dt, dposition, preliminary=.true.)
         call calc_force(position, force)
         if(abs(force)>0) then
            ! The code below cannot handle averaging spreading operators yet...
            stop "Spreading of force for particle_integrator==21 not implemented"
         end if   
      end if      

   end if

   ! Sample random numbers
   if(integrator>0) then ! We might need two RNGs per step
      call NormalRNGs(rngs, size(rngs))
   else      
      call NormalRNGs(rngs(:,1), size(rngs(:,1)))
   end if   

   old_density=density ! Save this
   select case(abs(integrator))
   case(0,10) ! Euler
      ! Euler step:
      !------------------------

      ! Explicit forcing:
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt)
         end if
      end if

      ! Stochastic forcing:
      massflux = rngs(:,1)
      call StochasticUpdate(density,massflux,ddensity)

      if(integrator==0) then ! Forward Euler
         call ExplicitUpdate(density,ddensity,weight=1.0_wp)
      else ! Backward Euler  
         call ImplicitUpdate(density,ddensity,weight=1.0_wp)
      end if   

      density(1:n_cells) = density(1:n_cells) + ddensity         
      call FixPeriodic(density)
      !------------------------

   case(1,11) ! Trapezoidal PC
      denscopy=density ! Make a copy of the old state

      ! Euler predictor at endpoint:
      !------------------------
      if(integrator>0) then ! 2RNG scheme is inconsistent
         massflux = sqrt(2.0_wp)*rngs(:,1)
      else ! Second-order
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            poscopy=position
            call move_particle(dt, dposition, preliminary=.true.)
         end if
      end if
      call StochasticUpdate(density,massflux,ddensity)
      ddensity_copy = ddensity ! Store explicit fluxes

      call DiffusiveUpdate(density,ddensity,1.0_wp)
      density(1:n_cells) = density(1:n_cells) + ddensity

      call FixPeriodic(density)
      !------------------------

      ! Trapezoidal corrector
      !------------------------
      if(integrator>0) then ! 2RNG scheme is inconsistent
         massflux = sqrt(2.0_wp)*rngs(:,2)
      else ! Second-order
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         if(particle_integrator==1) call calc_force(position, force)
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt,dposcopy, preliminary=.true.)
            position=poscopy+(dposition+dposcopy)/2 ! Trapezoidal
            call fix_particle()
         end if
      end if
      call StochasticUpdate(density,massflux,ddensity)

      if(abs(integrator)>=10) then ! Implicit diffusion        
         ddensity = 0.5_wp*(ddensity + ddensity_copy)
         call DiffusiveUpdate(denscopy,ddensity,1.0_wp)
         density(1:n_cells) = denscopy(1:n_cells) + ddensity            
      else ! Explicit
         call DiffusiveUpdate(density,ddensity,1.0_wp)
         density(1:n_cells) = 0.5_wp*(denscopy(1:n_cells) + &
                                      density(1:n_cells) + ddensity)
      end if

      call FixPeriodic(density)
      !------------------------

   case(2,12) ! Midpoint PC
      ! Note: We follow Boyce's algorithm: The predictor is Crank-Nicolson full time step
      denscopy=density ! Make a copy of the old state

      ! Euler predictor at midpoint:
      !------------------------
      if(integrator==2) then ! Explicit 2RNG scheme
         step_length=0.5_wp ! Step to the midpoint only
         massflux = sqrt(step_length)*rngs(:,1)
      else if(integrator==12) then ! 2RNG scheme should be second order   
         step_length=1.0_wp ! Step to the end in the predictor, then take midpoint
         massflux = sqrt(2.0_wp)*rngs(:,1) ! Should be weakly second-order
      else ! 1RNG scheme is not consistent
         step_length=1.0_wp ! Step to the end in the predictor, then take midpoint
         massflux = sqrt(step_length)*rngs(:,1)         
      end if

      call AdvectiveUpdate(density,ddensity,step_length)
      if(initcond<0) then
         call spread(force*step_length*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            poscopy=position
            call move_particle(dt/2, preliminary=.true.)
         end if            
      end if
      call StochasticUpdate(density,massflux,ddensity)
      call DiffusiveUpdate(density,ddensity,step_length)

      ! Estimate midpoint:
      if(integrator==2) then ! We are at the midpoint now
         density(1:n_cells) = density(1:n_cells) + ddensity
      else ! Average the beginning and end
         density(1:n_cells) = density(1:n_cells) + 0.5_wp*ddensity
      end if
      call FixPeriodic(density)
      !------------------------

      ! Midpoint corrector
      !------------------------
      if(integrator>0) then ! Add the Weiner increments from the two halves
         massflux = sqrt(0.5_wp)*rngs(:,1) + sqrt(0.5_wp)*rngs(:,2)
      else
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         if(particle_integrator==1) call calc_force(position, force)
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt, dposition, preliminary=.true.)
            position = poscopy + dposition
            call fix_particle()
         end if
      end if
      call StochasticUpdate(density,massflux,ddensity)

      if(abs(integrator)>=10) then ! Implicit diffusion        
         call DiffusiveUpdate(denscopy,ddensity,1.0_wp)
      else ! Explicit
         call DiffusiveUpdate(density,ddensity,1.0_wp)                        
      end if
      density(1:n_cells) = denscopy(1:n_cells) + ddensity
      call FixPeriodic(density)

      !------------------------

   case(13,14) ! L-stable trapezoidal (w2=1) or midpoint (w2=1/2)
      denscopy=density ! Make a copy of the old state

      ! Euler predictor at point w2*dt:
      !------------------------
      massflux = sqrt(w2)*rngs(:,1)
      call AdvectiveUpdate(density,ddensity,w2)
      if(initcond<0) then
         call spread(force*w2*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            poscopy=position
            call move_particle(w2*dt,dposition, preliminary=.true.)
            dposition=dposition/w2 ! Store F*dt for later
         end if
      end if
      ddensity_copy = ddensity/w2 ! Store advective fluxes for reuse
      call StochasticUpdate(density,massflux,ddensity)

      call ExplicitUpdate(density,ddensity,weight=(w2-w1)) ! Explicit diffusive terms
      call ImplicitUpdate(density,ddensity,weight=w1)

      density(1:n_cells) = density(1:n_cells) + ddensity         
      call FixPeriodic(density)
      !------------------------

      ! Trapezoidal corrector
      !------------------------
      if(integrator>0) then ! Add the Weiner increments from the two parts
         massflux = sqrt(w2)*rngs(:,1) + sqrt(1-w2)*rngs(:,2)
      else
         massflux = rngs(:,1)
      end if   
      ! Calculate new advective fluxes:
      call AdvectiveUpdate(density,ddensity,w5)
      if(initcond<0) then
         if(particle_integrator==1) call calc_force(position, force)
         call spread(force*w5*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit         
            call move_particle(w5*dt, dposcopy, preliminary=.true.)
            position = poscopy + dposcopy + (1-w5)*dposition
            call fix_particle()
         end if   
      end if
      ! Add the old ones with the right weights:
      ddensity = ddensity + (1-w5)*ddensity_copy

      call StochasticUpdate(density,massflux,ddensity)
      ! Explicit diffusive terms:
      call ExplicitUpdate(denscopy,ddensity,weight=(1-w3-w4))
      call ExplicitUpdate(density,ddensity,weight=w3)
      ! And now do the implicit solve:
      call ImplicitUpdate(denscopy,ddensity,weight=w4)

      density(1:n_cells) = denscopy(1:n_cells) + ddensity         
      call FixPeriodic(density)
      !------------------------

   case(15) ! Crank-Nicolson (diffusion) + Euler (rest) method 
      ! This is *not* second-order accurate!

      massflux = rngs(:,1)
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt)
         end if
      end if
      call StochasticUpdate(density,massflux,ddensity)

      call ExplicitUpdate(density,ddensity,weight=0.5_wp) ! Explicit diffusive terms
      call ImplicitUpdate(density,ddensity,weight=0.5_wp)

      density(1:n_cells) = density(1:n_cells) + ddensity         
      call FixPeriodic(density)

   case(3,4,5) ! Explicit RK3
      denscopy=density ! Make a copy of the old state

      ! RK stage 1:
      !------------------------
      if(integrator>0) then
         massflux = rngs(:,1) + w1 * rngs(:,2)            
      else
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            poscopy=position
            call move_particle(dt, dposition, preliminary=.true.)
         end if            
      end if
      call StochasticUpdate(density,massflux,ddensity)
      call DiffusiveUpdate(density,ddensity,1.0_wp)
      density(1:n_cells) = density(1:n_cells) + ddensity ! Simple Euler
      call FixPeriodic(density)

      ! RK stage 2:
      !------------------------
      if(integrator>0) then
         massflux = rngs(:,1) + w2*rngs(:,2)
      else
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         if(particle_integrator==1) call calc_force(position, force)
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt, dposcopy, preliminary=.true.)
            position = poscopy + (dposition+dposcopy)/4
            dposition=(dposition+dposcopy)/4 ! Store for stage 3
            call fix_particle()
         end if               
      end if
      call StochasticUpdate(density,massflux,ddensity)
      call DiffusiveUpdate(density,ddensity,1.0_wp)
      density(1:n_cells) = density(1:n_cells) + ddensity
      call FixPeriodic(density)

      density = (3*denscopy+density)/4 ! Stage 2 weighted average

      ! RK stage 3:
      !------------------------
      if(integrator>0) then
         massflux = rngs(:,1) + w3 * rngs(:,2)            
      else
         massflux = rngs(:,1)
      end if   
      call AdvectiveUpdate(density,ddensity,1.0_wp)
      if(initcond<0) then
         if(particle_integrator==1) call calc_force(position, force)
         call spread(force*dt, ddensity, position)
         if(particle_integrator==1) then ! Explicit
            call move_particle(dt, dposcopy, preliminary=.true.)
            position = poscopy + 2*(dposition+dposcopy)/3
            call fix_particle()
         end if               
      end if
      call StochasticUpdate(density,massflux,ddensity)
      call DiffusiveUpdate(density,ddensity,1.0_wp)
      density(1:n_cells) = density(1:n_cells) + ddensity
      call FixPeriodic(density)

      density = (denscopy+2*density)/3 ! Stage 3 weighted average
      !------------------------

   case default
      stop "Time stepping scheme not recognized!"
   end select   

   if(initcond<0) then ! Finish the particle move
            
      select case(particle_integrator)
      case(0) ! Use Lie splitting         
         call move_particle(dt, dposition)
         !write(41,*) time, dposition**2/(2*dt), dposition
      case(-1,-2,-3) ! Use Strang splitting
         call move_particle(dt/2)
      case(21) ! Boyce's trapezoidal scheme
         call move_particle(dt, dposcopy, preliminary=.true.)
         position=poscopy+(dposition+dposcopy)/2
         call fix_particle()      
      case(22) ! Boyce's midpoint scheme
         denscopy=density
         density=(density+old_density)/2 ! Midpoint estimate
         call move_particle(dt, dposition, preliminary=.true.)
         position = poscopy + dposition
         call fix_particle()
         density=denscopy
      case(10,11,12) ! Staggered Euler/Midpoint
         ! We are moving from step n to n+1, so now
         ! pos_stag_last=q(n+1/2), pos_stag_next=q(n+3/2)
      
         ! We need to move the particle to (n+3/2)
         position_shift=periodic_shift
         call move_particle(dt, dposition)

         pos_stag_last=pos_stag_next ! n+1/2
         pos_stag_next=position ! n+3/2

         ! Estimate non-staggered position at (n+1) as midpoint:
         position = pos_stag_last + dposition/2
         
         shifted_position = position + position_shift ! Remember this before continuing
         call fix_particle(preliminary=.true.)
         !write(32,*) time, real(pos_stag_last+position_shift), real(pos_stag_next+periodic_shift)

      end select

      select case(particle_integrator)
      case(10,11,12) ! Staggered Euler/Midpoint
         ! The shifting has to be done carefully here
      case default
         position_shift=periodic_shift
         shifted_position=position+position_shift
      end select

   end if
                        
end subroutine

subroutine InitializeFluid()
   
   if(periodic) then
      system_length=(n_cells+1)*dx
   else
      system_length=n_cells*dx
   end if   
   
   select case(initcond)
   case(0)
      ! Initialize with constant profile
      do cell=0, n_cells+1
         density(cell)=bcs(1) ! No gradient
      end do
   case(1)
      ! Initialize with linear profile
      do cell=0, n_cells+1
         density(cell)=bcs(1) + (bcs(2)-bcs(1))*real(cell,wp)/(n_cells+1) ! Gradient
      end do
   case(2,3)
      
      if(initcond==2) then ! Moving shock starts 1/4 into the domain
         shockpos=0.25_wp
      else ! Stationary shock starts in middle
         shockpos=0.5_wp
      end if  
            
      do cell=0, n_cells+1
         r_tmp=(cell*dx - shockpos*((n_cells+1)*dx)) ! x-x0
         density(cell) = 0.5_wp*(bcs(1)+bcs(2)) - 0.5_wp*(bcs(1)-bcs(2))* &
            tanh(r_tmp*abs(advspeed)*(bcs(1)-bcs(2))/(4*diffcoeff))            
      end do
      ! call write_shock()

   case(4) ! Periodic "triangle-wave" like solution
      
      if(.not.periodic) then
         stop "Must have periodic BCs for initcond==4"
      end if       
      
      !time=0; call periodic_solution_simple(density,time) ! This works for any size domain
      ! This requires length=2*Pi
      time=0.25_wp; call periodic_solution(density,time) ! Starts with nearly triangular wave
      call write_solution(density)
   
   case default
   
      ! Initialize at equilibrium
      if(.false.) then
         do cell=0, n_cells+1
            density(cell)=bcs(1) ! No gradient
         end do
      else
         do cell=0, n_cells+1
            density(cell)=bcs(1) + (bcs(2)-bcs(1))*real(cell,wp)/(n_cells+1) ! Gradient
            !density(cell)=bcs(1) + (bcs(2)-bcs(1))*(real(cell,wp)/(n_cells+1))**2 ! Quadratic
         end do      
      end if
         
   end select
   call FixPeriodic(density)

end subroutine

subroutine FixPeriodic(density)
   real(wp), intent(inout) :: density(0:n_cells+1)
   if(periodic) then
      density(n_cells+1)=density(1)
      density(0)=density(n_cells)
   end if
end subroutine

function A_u(u) ! Amplitude of stochastic flux
   real(wp), intent(in) :: u
   real(wp) :: A_u

   select case(abs(noise))
   case(0)
      A_u=1.0_wp
   case(1)
      A_u=u
   case(2)
      A_u=(1-u)/2
   case(3)
      A_u=u*(1-u)
   case(4)
      A_u=(1-u**2)/4
   case(5)
      A_u=u**2
   case default
      A_u=0   
   end select   

end function

subroutine AdvectiveUpdate(density,ddensity,step_length) ! Calculate A(u)
   real(wp), intent(in) :: density(0:n_cells+1)
   real(wp), intent(out) :: ddensity(1:n_cells) ! The explicit flux update
   real(wp), intent(in) :: step_length ! 0.5 for midpoint, 1 otherwise
   
   ! u_t + v*u_x = D*u_xx + [epsilon*2*D*A(u)*W]_x
   ! v=c*u if c>0 (Burgers), or v=abs(c) if c<=0 (advection-diffusion)
   do cell=1, n_cells
   
      ddensity(cell)=0.0_wp

      ! Calculate advective fluxes by first calculating f_right-f_left
      if(advspeed<=0) then ! Advection-diffusion equation (linear flux=u)
         ! Simple centered difference for u_x:
         ddensity(cell) = 0.5_wp*(density(cell+1)-density(cell-1))
      else ! Burgers equation (nonlinear flux=u^2/2)
      
         select case(interpolations(1))
         case(0) ! Simple average
            ! u_(j+1/2)=(u_j + u_(j+1))/2
            fdensity = 0.5_wp*(density(cell+1)+density(cell)) ! Right cell-interpolated density
            ddensity(cell) = ddensity(cell) + (0.5_wp*fdensity**2)
            fdensity = 0.5_wp*(density(cell-1)+density(cell)) ! Left cell-interpolated density
            ddensity(cell) = ddensity(cell) - (0.5_wp*fdensity**2)               
         case(1) ! Hamiltonian-structure-preserving scheme
            ! u_(j+1/2)^2 = (u_j^2 + u_j*u_(j+1) + u_(j+1)^2)/3
            fdensity = (density(cell+1)**2 + density(cell+1)*density(cell) + density(cell)**2)/3
            ddensity(cell) = ddensity(cell) + (0.5_wp*fdensity)
            fdensity = (density(cell-1)**2 + density(cell-1)*density(cell) + density(cell)**2)/3
            ddensity(cell) = ddensity(cell) - (0.5_wp*fdensity)
         end select
               
      end if      
      ddensity(cell) = - step_length*dt*abs(advspeed)*ddensity(cell)/dx

      if(initcond==3) then ! Add advection term to make shock frozen in place
         ddensity(cell) = ddensity(cell) + step_length*dt*abs(advspeed)/dx* &
            0.25_wp*(bcs(1)+bcs(2))*(density(cell+1)-density(cell-1))
      end if   

   end do
   
end subroutine

! For explicit
! ddensity = (I+L*dt)*density + ddensity
! For implicit:
! ddensity = (I-L*dt/2)^(-1) * [ (I+L*dt/2)*density + ddensity ]
subroutine DiffusiveUpdate(density,ddensity,step_length) ! Linear (Laplacian) terms
   real(wp), intent(in) :: density(0:n_cells+1) ! Current state
   real(wp), intent(inout) :: ddensity(1:n_cells) ! The explicit part of the update (rhs)
   real(wp), intent(in) :: step_length ! Implicit system time step is step_length*dt
   
   ! Local variables for semi-implicit solve
   real(wp), dimension(1:n_cells) :: a,b,c,r
   real(wp) :: u(0:n_cells+1) ! Temporary copy of density for implicit solve
   
   do cell=1, n_cells ! Calculate diffusive fluxes explicitly
      r(cell) = ddensity(cell) + step_length*dt* &
         diffcoeff/dx**2 * (density(cell-1)-2*density(cell)+density(cell+1))                 
   end do
   
   if(abs(integrator)<10) then ! Explicit diffusion    
   
      ddensity = r   
   
   else ! Implicit midpoint (Crank-Nicolson) diffusion
   
      ! Off-diagonal entries
      a = -0.5_wp*step_length*dt*diffcoeff/dx**2
      ! Diagonal
      b = 1.0_wp + step_length*dt*diffcoeff/dx**2

      ! Now solve the tridiagonal system for u=density(t+step_length*dt)-density(t)
      call SolveTridiagonal(a,b,a, r, u(1:n_cells), n_cells, merge(1,0,periodic))
      if(periodic) then
         u(n_cells+1)=u(1)
         u(0)=u(n_cells)
      else ! Dirichlet BCs
         u(0)=0
         u(n_cells+1)=0
      end if
      
      ! Check the linear solver:
      if(.false.) then
         do cell=1, n_cells
            c(cell) = u(cell) - 0.5_wp*step_length*dt*&
                      diffcoeff/dx**2 * (u(cell-1)-2*u(cell)+u(cell+1))
         end do
         write(*,*) "MAX ERROR=", maxval(abs(c-r(1:n_cells)))
      end if
      
      ddensity = u(1:n_cells)

   end if

end subroutine

! For general IMEX-schemes we want more flexibility, so split into two steps:
! ddensity = ddensity+w*L*dt*density
subroutine ExplicitUpdate(density,ddensity,weight) ! Linear (Laplacian) terms
   real(wp), intent(in) :: density(0:n_cells+1) ! Current state
   real(wp), intent(inout) :: ddensity(1:n_cells) ! The explicit fluxes
   real(wp), intent(in) :: weight ! Time step is weight*dt
      
   do cell=1, n_cells ! Calculate diffusive fluxes explicitly
      ddensity(cell) = ddensity(cell) + weight*dt* &
         diffcoeff/dx**2 * (density(cell-1)-2*density(cell)+density(cell+1))                 
   end do

end subroutine

! density = (I-weight*L*dt)^(-1) * [ density + ddensity ]
! We rewrite this in terms of increments as:
! ddensity = (I-weight*L*dt)^(-1) * [ weight*L*dt*density + ddensity ]
subroutine ImplicitUpdate(density,ddensity,weight) ! Linear (Laplacian) terms
   real(wp), intent(in) :: density(0:n_cells+1) ! Current state
   real(wp), intent(inout) :: ddensity(1:n_cells) ! The explicit part of the update (rhs)
   real(wp), intent(in) :: weight ! Implicit system time step is weight*dt
   
   ! Local variables for semi-implicit solve
   real(wp), dimension(1:n_cells) :: a,b,c,r
   real(wp) :: u(0:n_cells+1) ! Temporary copy of density for implicit solve
   
   do cell=1, n_cells ! Calculate rhs term weight*L*dt*density explicitly
      r(cell) = ddensity(cell) + weight*dt* &
         diffcoeff/dx**2 * (density(cell-1)-2*density(cell)+density(cell+1))            
   end do
   
   ! Off-diagonal entries
   a = -weight*dt*diffcoeff/dx**2
   ! Diagonal
   b = 1.0_wp + 2*weight*dt*diffcoeff/dx**2

   ! Now solve the tridiagonal system for u=density(t+weight*dt)-density(t)
   call SolveTridiagonal(a,b,a, r, u(1:n_cells), n_cells, merge(1,0,periodic))
   if(periodic) then
      u(n_cells+1)=u(1)
      u(0)=u(n_cells)
   else ! Dirichlet BCs
      u(0)=0
      u(n_cells+1)=0
   end if

   ! Check the linear solver:
   if(.false.) then
      do cell=1, n_cells
         c(cell) = u(cell) - weight*dt*&
                   diffcoeff/dx**2 * (u(cell-1)-2*u(cell)+u(cell+1))
      end do
      write(*,*) "MAX ERROR=", maxval(abs(c-r(1:n_cells)))
   end if

   ddensity = u(1:n_cells)

end subroutine

subroutine StochasticUpdate(density,massflux,ddensity)
   real(wp), intent(in) :: density(0:n_cells+1)
   real(wp), intent(inout) :: massflux(1:n_cells+1) ! Comes in filled with random numbers
   real(wp), intent(inout) :: ddensity(1:n_cells) ! Add divergence of the stochastic fluxes to this

   if(noise>0) then ! Multiplicative noise
   do face=1, n_cells+1
      select case(interpolations(2))
      case(0)
        fdensity = 0.5_wp*(density(face-1)+density(face)) ! Face-interpolated density
        stochforcing(face) = sqrt(2*A_u(fdensity))
      case(1)  
        stochforcing(face) = sqrt((A_u(density(face-1))+A_u(density(face))))
      case(2)
         if(abs(interpolations(2))/=3) stop "Fancy interpolation only implemented for noise=3"
         stochforcing(face) = sqrt((density(face-1)*(1-density(face)) + &
                                    density(face)*(1-density(face-1))))
      case(3) ! Split fluxes into left-to-right and right-to-left                                                 
         if(abs(noise)/=3) stop "Directional interpolation only implemented for noise=3"
         call NormalRNG(l2r) ! Left to right flux
         call NormalRNG(r2l) ! Right to left flux
         massflux(face) = sqrt(density(face-1)*(1-density(face)))*l2r - &
                          sqrt(density(face)*(1-density(face-1)))*r2l
      end select
   end do
   end if

   if((noise>0).and.(interpolations(2)==3)) then
      massflux = massflux*sqrt(varsc)
   else
      massflux = massflux*stochforcing*sqrt(varsc)      
   end if   
   if(periodic) then
      massflux(n_cells+1)=massflux(1)
   end if

   do cell=1, n_cells
      ! And stochastic fluxes:
      ddensity(cell) = ddensity(cell) + sqrt(diffcoeff*abs(dt)/dx**3) * (massflux(cell+1)-massflux(cell))
   end do

end subroutine

!======================================================
! Immersed-particle tools
!======================================================
subroutine move_particle(dtime, dposition, preliminary)
   real(wp), intent(in) :: dtime
   real(wp), intent(out), optional :: dposition
   logical, optional, intent(in) :: preliminary

   real(wp) :: position_old
   
   call interpolate(velocity, density, position)

   !write(102,*) time, real(velocity)

   position_old=position
   select case(particle_integrator)
   case(-1,11) ! Trapezoidal
      position = position + dtime*velocity
      call fix_particle(preliminary=.true.)
      call interpolate(velocity_next, density, position)
      position = position_old + dtime*(velocity+velocity_next)/2
   case(-2,12) ! Midpoint
      position = position + dtime/2*velocity
      call fix_particle(preliminary=.true.)
      call interpolate(velocity_next, density, position)
      position = position_old + dtime*velocity_next
   case default ! Euler
      position = position + dtime*velocity
   end select
   ! We need to calculate this before doing periodic BCs:
   if(present(dposition)) dposition=position-position_old
   
   call fix_particle(preliminary)

end subroutine

subroutine fix_particle(preliminary) ! Fix periodic BCs for particle   
   logical, optional, intent(in) :: preliminary
   
   logical :: preliminary_
   
   preliminary_=.false.
   if(present(preliminary)) preliminary_=preliminary

   ! Make sure the particle stays inside the valid domain
   if(periodic) then
      if(position<0) then
         position=position+n_cells*dx
         if(.not.preliminary_) then
            periodic_shift=periodic_shift-n_cells*dx
         end if
      else if(position>n_cells*dx) then
         position=position-n_cells*dx
         if(.not.preliminary_) then
            periodic_shift=periodic_shift+n_cells*dx
         end if   
      end if
   else if(position<dx) then
      stop "Particle went too far to the left"
   else if(position>n_cells*dx) then
      stop "Particle went too far to the right"
   end if   
end subroutine   

subroutine calc_force(position, force)
   real(wp), intent(in) :: position ! q
   real(wp), intent(out) :: force ! F(q)
   
   select case(initcond)
   case(-2)
      force=1.0_wp
   case(-3)
      force=-(position-0.5_wp*n_cells*dx)/(0.5_wp*n_cells*dx)
   case default
      force=0.0_wp   
   end select
   force = force_const*force
   
end subroutine

subroutine interpolate(velocity, density, position)
   real(wp), intent(in) :: density(0:n_cells+1) ! v
   real(wp), intent(in) :: position ! Particle position q inside [0,Length]
   real(wp), intent(out) :: velocity ! u=J(q)*v

   integer :: cells(4)
   real(wp) :: weights(4)
   
   ! We can rely here on the ghost cells but there is only one
   cells(2)= floor(position/dx)
   cells(1)=cells(2)-1
   cells(3)=cells(2)+1
   cells(4)=cells(2)+2

   ! On the left we can go outside of the system
   if(periodic) then
      if(cells(1)<0) cells(1)=cells(1)+n_cells
      if(cells(4)>n_cells+1) cells(4)=cells(4)-n_cells
   end if   
      
   if(any(cells<0).or.any(cells>n_cells+1)) then
      write(*,*) "Bad position=", position
      stop "Trying to intepolate outside of valid domain"
   end if
   
   call approx_dirac(distance=(cells(3)*dx-position)/dx, npt=kernel_width, w=weights)
   
   velocity=sum(weights*density(cells))

   if(.false.) then
      write(*,*) "position=", position
      write(*,*) "cells=", cells
      write(*,*) "weights=", weights
      write(*,*) "velocities=", real(cells)/(n_cells+1)
      write(*,*) "velocities=", real(cells)/(n_cells+1)
      write(*,*) real(velocity), real(position/((n_cells+1)*dx))
   end if

end subroutine

! The operator S is not dimensionless: It has units 1/dx:
subroutine spread(force, fdensity, position) ! No ghost cells in fdensity!
   real(wp), intent(inout) :: fdensity(1:n_cells) ! f=f+S(q)*F
   real(wp), intent(in) :: position ! Particle position q
   real(wp), intent(in) :: force ! Force F

   integer :: cells(4)
   real(wp) :: weights(4)
   
   cells(2)= floor(position/dx)
   cells(1)=cells(2)-1 
   cells(3)=cells(2)+1
   cells(4)=cells(2)+2
   
   ! We do not spread force to the ghost cells
   if(periodic) then ! Spread the force to the interior
      if(cells(1)<1) cells(1)=cells(1)+n_cells
      if(cells(2)<1) cells(2)=cells(2)+n_cells
      if(cells(4)>n_cells) cells(4)=cells(4)-n_cells
   else ! Abort if too close to the boundary
      if(cells(1)<1) stop "Spreading force to boundary cells < 0"
      if(cells(4)>n_cells) stop "Spreading force to boundary cells > 0"
   end if
   
   call approx_dirac(distance=(cells(3)*dx-position)/dx, npt=kernel_width, w=weights)
   
   fdensity(cells)=fdensity(cells)+force*weights/dx
   if(periodic.and.(initcond==-2)) then
      ! Keep mean force zero to keep the system from accelerating
      fdensity(1:n_cells)=fdensity(1:n_cells)-force/dx/n_cells
   end if
   
end subroutine

subroutine approx_dirac(distance, npt, w)
   real(wp), intent(in) :: distance ! Must be in [0-1]
   integer, intent(in) :: npt
   real(wp), intent(out) :: w(4)
   
   real(wp) :: r, q
   
   r=abs(distance)
   
   if ((npt==4) .and. (r<=1.0_wp)) then ! Use the 4pt function: 0<=q<1   
      q = sqrt(1 + 4*r*(1-r))
      w(1) = (1 + 2*r - q)/8 ! At r-2
      w(2) = (1 + 2*r + q)/8 ! At r-1
      w(3) = (3 - 2*r + q)/8 ! At r
      w(4) = (3 - 2*r - q)/8 ! At r+1
   elseif ((npt==3) .and. (r<=0.5_wp)) then ! Use the 3pt function: 0<=r<1/2
      q = sqrt(1 - 3*r**2)
      w(1)=0
      w(2) = (2 + 3*r - q)/6 ! At r-1
      w(3) = (1 + q)/3 ! At r
      w(4) = (2 - 3*r - q)/6 ! At r +1
   elseif ((npt==3) .and. (r<=1.0_wp)) then ! Use the 3pt function: 1/2<=r<1
      r=1-r
      q = sqrt(1 - 3*r**2)
      w(4) = 0
      w(3) = (2 + 3*r - q)/6
      w(2) = (1 + q)/3
      w(1) = (2 - 3*r - q)/6
   elseif ((npt==2) .and. (r<=1.0_wp)) then ! Use the 2pt function: 0<=r<1
      w(1) = 0 ! At r-2
      w(2) = r ! At r-1
      w(3) = 1-r ! At r
      w(4) = 0   
   else
      !call abort()
      stop "Problem in approx_dirac"
   end if
   
end subroutine

!======================================================
! Exact solutions and tools for Burgers
!======================================================

! This is obtained from the Cole-Hopf transform
subroutine periodic_solution_simple(density, time)
   real(wp), intent(inout) :: density(0:n_cells+1)
   real(wp), intent(in) :: time ! Current time
   
   real(wp) :: x, n, a, b
   
   a=1.0_wp
   b=0.9_wp
   
   if(advspeed<0) stop "Advective speed must be positive for periodic_solution"

   n=2*pi/(n_cells*dx)
   do cell=1, n_cells
      x=cell*dx ! Goes from 0 to 2*Pi
      density(cell) = 2*diffcoeff*b*n*sin(n*x)/advspeed/ &
         (a*exp(diffcoeff*n**2*time)+b*cos(n*x))
   end do
   density(0)=density(n_cells)
   density(n_cells+1)=density(1)
   
end subroutine

! This is an infinite series solution by Fay that starts as a triangle wave
! Note that this one only works for (-Pi,Pi)
subroutine periodic_solution(density, time)
   real(wp), intent(inout) :: density(0:n_cells+1)
   real(wp), intent(in) :: time ! Current time
   
   real(wp) :: x, factor
   integer :: term

   if(advspeed<0) stop "Advective speed must be positive for periodic_solution"

   do cell=1, n_cells
      x=(cell-0.5_wp*n_cells)*dx ! In (-Pi,Pi)
      density(cell)=0
      do term=1, 100 ! One hundred terms should be enough
         ! Nonlinear equation:
         density(cell) = density(cell) - 2*diffcoeff/advspeed*&
            sin(term*x)/sinh(term*diffcoeff*time)
         ! Linear equation (note the n^2) in the denominator
         !density(cell) = density(cell) - sin(term*x)/exp(term**2*diffcoeff*time)
      end do   
   end do
   density(0)=density(n_cells)
   density(n_cells+1)=density(1)
   
end subroutine

subroutine write_solution(exact)
   real(wp), intent(in) :: exact(0:n_cells+1)

   write(15,*) "# t=", time
   write(16,*) "# t=", time
   do cell=merge(1,0,periodic), merge(n_cells,n_cells+1,periodic)
      write(15,*) cell*dx, density(cell), exact(cell)
      write(16,*) cell*dx, density(cell) - exact(cell)
   end do
   write(15,*) ! Blank line
   write(16,*) ! Blank line

end subroutine

subroutine write_shock()

   write(15,*) "# t=", time
   if(initcond==2) write(16,*) "# t=", time
   do cell=0, n_cells+1
      if(initcond==2) then ! Moving shock
         r_tmp = ( cell*dx - 0.25_wp*((n_cells+1)*dx)- &
            0.5_wp*abs(advspeed)*(bcs(1)+bcs(2))*time ) ! x-x0-s*t
      else
         r_tmp = cell*dx - 0.5_wp*((n_cells+1)*dx)
      end if      
      r_tmp = 0.5_wp*(bcs(1)+bcs(2)) - 0.5_wp*(bcs(1)-bcs(2))* &
         tanh(r_tmp*abs(advspeed)*(bcs(1)-bcs(2))/(4*diffcoeff))
      write(15,*) cell*dx, density(cell), r_tmp
      if(initcond==2) write(16,*) cell*dx, density(cell) - r_tmp
   end do
   write(15,*) ! Blank line
   if(initcond==2) write(16,*) ! Blank line

end subroutine

subroutine fit_shock(mincenter, shockpos)
   real(wp), intent(out) :: mincenter, shockpos

   real(wp) :: ubar, center, r_tmp, norm, minnorm
   integer :: shift, width
   
   ubar=sum(density(1:n_cells))/n_cells   
   shockpos = (ubar - bcs(2))/(bcs(1)-bcs(2))*n_cells*dx
   
   ! Search a small neighbourhood around the center to find best position
   minnorm=huge(1.0_wp)
   width=min(n_cells/4, max(16,n_cells/16))
   !write(*,*) "Searching width:", width
   do shift=-width, width
      if(.true.) then ! To avoid bias, add random shift
         call UniformRNG(r_tmp)
         r_tmp=r_tmp-0.5_wp
         center = shockpos + (shift + r_tmp)*dx
      else   
         center = shockpos + 0.25_wp*shift*dx ! Search to within accuracy dx/2
      end if   

      norm=0.0_wp ! L2 norm
      do cell=1, n_cells
         r_tmp = cell*dx - center
         r_tmp = 0.5_wp*(bcs(1)+bcs(2)) - 0.5_wp*(bcs(1)-bcs(2))* &
            tanh(r_tmp*abs(advspeed)*(bcs(1)-bcs(2))/(4*diffcoeff))
         norm = norm + (density(cell)-r_tmp)**2         
      end do
      !write(*,*) shift, norm, center
      
      if(norm<minnorm) then
         minnorm=norm
         mincenter=center
      end if

   end do
   
   !write(*,*) "Estimated shock positions: ", real(shockpos), real(mincenter)
   !write(*,*) "Shift in shock position = ", real((mincenter-shockpos)/dx)
   
end subroutine

subroutine write_fitted_shock(mincenter, shockpos)
   real(wp), intent(in) :: mincenter, shockpos

   real(wp) :: center, r_tmp, fit(2)
   integer :: shift, width
   
   shift=0
   do cell=1, n_cells
      r_tmp = cell*dx - (shockpos + 0.25_wp*shift*dx)      
      fit(1) = 0.5_wp*(bcs(1)+bcs(2)) - 0.5_wp*(bcs(1)-bcs(2))* &
         tanh(r_tmp*abs(advspeed)*(bcs(1)-bcs(2))/(4*diffcoeff))

      r_tmp = cell*dx - (mincenter + 0.25_wp*shift*dx)      
      fit(2) = 0.5_wp*(bcs(1)+bcs(2)) - 0.5_wp*(bcs(1)-bcs(2))* &
         tanh(r_tmp*abs(advspeed)*(bcs(1)-bcs(2))/(4*diffcoeff))

      write(17,*) cell*real(dx), real(density(cell)), real(fit)   
   end do
   write(17,*)

end subroutine

end program
