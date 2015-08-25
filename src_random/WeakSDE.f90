module IntegrateSDEs
   use Precision
   use BoxLibRNGs
   implicit none
   
   integer, parameter :: wp=r_dp, ip=i_sp ! Working precision

contains

   function f(x,t)
      real(wp) :: x,t,f
      f=x**2-x**3
   end function   

   function g(y)
      real(wp) :: y,g
      g=y+y**2-y**3
   end function

   subroutine ExplicitIntegrator(x,t0,x0,epsilon,dt,n_steps,method)
      ! Solves the one-dimension SDE:
      ! x'=f(x,t)+epsilon*W(t)
      real(wp), intent(out) :: x
      real(wp), intent(in) :: t0,x0,epsilon,dt
      integer, intent(in) :: n_steps,method
      
      integer :: step
      real(wp) :: t, RNG1, RNG2, w1, w2, w3
      real(wp) :: f_start, x_mid, x_end, f_end, x1, x2, x3
      
      x=x0
      t=t0
      do step=1, n_steps
         call NormalRNG(RNG1)
         
         select case(abs(method))
         case(0) ! Euler
            x = x + f(x,t)*dt + epsilon*sqrt(dt)*RNG1
            t = t + dt
         case (1) ! Trapezoidal
            f_start=f(x,t)
            x_end = x + f_start*dt + epsilon*sqrt(dt)*RNG1 ! Euler predictor
            t = t + dt
            if(method>0) then ! 2RNGs (first-order)
               call NormalRNG(RNG2)
               x = x + (f_start+f(x_end,t))/2*dt + epsilon*sqrt(dt)*RNG2
            else ! 1RNG (second-order)
               x = x + (f_start+f(x_end,t))/2*dt + epsilon*sqrt(dt)*RNG1
            end if
         case(2) ! Midpoint
            x_mid = x + f(x,t)*dt/2 + epsilon*sqrt(dt/2)*RNG1
            t = t + dt/2
            if(method>0) then ! 2RNGs (second-order)
               call NormalRNG(RNG2)
               x = x + f(x_mid,t)*dt + epsilon*sqrt(dt/2)*(RNG1+RNG2)
            else ! 1RNG (first-order)
               x = x + f(x_mid,t)*dt + epsilon*sqrt(dt)*RNG1
            end if
            t = t + dt/2
         case(3,4) ! RK3
         
            call NormalRNG(RNG2)
            
            if(method==4) then ! Old bad weights
               w1=-sqrt(3.0_wp)
               w2=sqrt(3.0_wp)
               w3=0.0_wp
            else if(method>0) then ! 2RNGs (second-order)
               ! New 2RNG scheme:
               w1=(2*sqrt(2.0_wp)-sqrt(3.0_wp))/5
               w2=(-4*sqrt(2.0_wp)-3*sqrt(3.0_wp))/5
               w3=(sqrt(2.0_wp)+2*sqrt(3.0_wp))/10
            else ! 1RNG (first-order)
               w1=0; w2=0; w3=0;
            end if   
            
            x1 = x + f(x,t)*dt + epsilon*sqrt(dt)*(RNG1+w1*RNG2)
            x2 = 3*x/4 + (x1 + f(x1,t+dt)*dt + epsilon*sqrt(dt)*(RNG1+w2*RNG2))/4
            x3 = x/3 + 2*(x2 + f(x2,t+dt/2)*dt + epsilon*sqrt(dt)*(RNG1+w3*RNG2))/3
            
            x=x3      
            t = t + dt            
         case default
            stop "Unregnized SDE method"                              
         end select   
      end do
   end subroutine

end module

program WeakSDE
   use Precision
   use IntegrateSDEs
   use Simple_Sorting
   use iso_c_binding
   
   interface GSLSort ! GNU Scientific Library heapsort
      !void gsl_sort (double * data, size_t stride, size_t n)
      subroutine gsl_sort(data, stride, size) bind(c)
         use iso_c_binding
         integer(c_size_t), value :: stride, size
         real(c_double), intent(inout) :: data(size)
      end subroutine
      subroutine gsl_sort_float(data, stride, size) bind(c)
         ! Uses heapsort
         use iso_c_binding
         integer(c_size_t), value :: stride, size
         real(c_float), intent(inout) :: data(size)
      end subroutine      
   end interface
   
   ! n_runs=10000000; L1_std=0.0003;   L2_std=0.0015;   L12_noise=0.0004;
   integer(i_dp) :: run
   integer(i_dp), parameter :: n_runs=1600000 !00
   real(wp) :: x0, T, epsilon, x, t0, dt
   real(wp) :: dts, errors_L1, errors_L2
   integer :: min_pow, max_pow, dt_pow, method
   
   !real(wp) :: results(n_runs), exact_distro(n_runs)
   ! It is important to save memory here:
   real(r_sp), dimension(:), allocatable :: results, exact_distro
   
   allocate(results(n_runs), exact_distro(n_runs))

   x0=0
   t0=0
   T=1 ! Time of approximation
   epsilon=1.0   

   min_pow=3
   max_pow=8

   read(*,*) method

   Change_dt: do dt_pow=max_pow,min_pow,-1
      dt=1.0_wp/2**dt_pow
      !dt=1.0_wp/2**4 ! Temporary testing
      n_steps=nint(T/dt)
      !write(*,*) "Running ", n_steps, "x", n_runs, " steps"

      do run=1,n_runs
         call ExplicitIntegrator(x,t0,x0,epsilon,dt,n_steps,method)
         results(run)=g(x)         
      end do
      write(*,*) "<g(x)>=", sum(results)/n_runs
      
      if(dt_pow==max_pow) then ! Most-accurate run
         exact_distro=results
         call HeapSort(exact_distro) ! Too slow and recursive
         !call GSLSort(exact_distro, 1_c_size_t, int(n_runs,c_size_t))
         !write(*,*) "Sorted=", exact_distro
         cycle Change_dt
      end if

      ! Compute L1 and L2 Wesserstein distance between distros
      results=results
      call HeapSort(results)
      !call GSLSort(results, 1_c_size_t, int(n_runs,c_size_t))
      errors_L1 = (sum(abs(results-exact_distro))/n_runs)
      errors_L2 = sqrt(sum((results-exact_distro)**2)/n_runs)
      
      write(*,*) "Errors=", real((/dt, errors_L1, errors_L2/))
      write(10+method,*) real((/dt, errors_L1, errors_L2/))

   end do Change_dt

end program
