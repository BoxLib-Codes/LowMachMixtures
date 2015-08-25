MODULE FFTW
   USE Precision, ONLY : r_dp, r_sp, i_sp, i_dp
   PUBLIC
   
   PRIVATE :: r_dp, r_sp, i_sp, i_dp ! Do not export these
   
   INTEGER, PARAMETER :: fftw_p=i_dp ! Plan kind
   
   INTERFACE FFTW_PlanDFT
      SUBROUTINE dfftw_plan_dft_1d(plan,n,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, sign, flags
         COMPLEX(KIND=r_dp), DIMENSION(n), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(n), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_2d(plan,nx,ny,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, sign, flags
         COMPLEX(KIND=r_dp), DIMENSION(nx,ny), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(nx,ny), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_3d(plan,nx,ny,nz,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, sign, flags
         COMPLEX(KIND=r_dp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(nx,ny,nz), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      
      SUBROUTINE sfftw_plan_dft_1d(plan,n,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, sign, flags
         COMPLEX(KIND=r_sp), DIMENSION(n), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(n), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_2d(plan,nx,ny,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, sign, flags
         COMPLEX(KIND=r_sp), DIMENSION(nx,ny), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(nx,ny), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_3d(plan,nx,ny,nz,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, sign, flags
         COMPLEX(KIND=r_sp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(nx,ny,nz), INTENT(OUT), TARGET :: out
      END SUBROUTINE
   END INTERFACE
   
   INTERFACE FFTW_PlanFT   
      SUBROUTINE dfftw_plan_dft(plan,rank,n,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, sign, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         COMPLEX(KIND=r_dp), DIMENSION(PRODUCT(n)), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: out
      END SUBROUTINE   
      SUBROUTINE sfftw_plan_dft(plan,rank,n,in,out,sign,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, sign, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         COMPLEX(KIND=r_sp), DIMENSION(PRODUCT(n)), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: out
      END SUBROUTINE   
   END INTERFACE

   INTERFACE FFTW_PlanDR2C ! FFT of real data (only half of the complex plane is stored)
      SUBROUTINE dfftw_plan_dft_r2c_1d(plan,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, flags
         REAL(KIND=r_dp), DIMENSION(n), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(n/2+1), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_r2c_2d(plan,nx,ny,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, flags
         REAL(KIND=r_dp), DIMENSION(nx,ny), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(nx/2+1,ny), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_r2c_3d(plan,nx,ny,nz,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, flags
         REAL(KIND=r_dp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION(nx/2+1,ny,nz), INTENT(OUT), TARGET :: out
      END SUBROUTINE

      SUBROUTINE sfftw_plan_dft_r2c_1d(plan,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, flags
         REAL(KIND=r_sp), DIMENSION(n), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(n/2+1), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_r2c_2d(plan,nx,ny,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, flags
         REAL(KIND=r_sp), DIMENSION(nx,ny), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(nx/2+1,ny), INTENT(OUT), TARGET :: out
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_r2c_3d(plan,nx,ny,nz,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, flags
         REAL(KIND=r_sp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION(nx/2+1,ny,nz), INTENT(OUT), TARGET :: out
      END SUBROUTINE
   END INTERFACE
      
   INTERFACE FFTW_PlanDC2R ! IFFT for real data (only half of the complex plane is stored)
      SUBROUTINE dfftw_plan_dft_c2r_1d(plan,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, flags
         REAL(KIND=r_dp), DIMENSION(n), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_dp), DIMENSION(n/2+1), INTENT(OUT), TARGET :: in
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_c2r_2d(plan,nx,ny,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, flags
         REAL(KIND=r_dp), DIMENSION(nx,ny), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_dp), DIMENSION(nx/2+1,ny), INTENT(OUT), TARGET :: in
      END SUBROUTINE
      SUBROUTINE dfftw_plan_dft_c2r_3d(plan,nx,ny,nz,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, flags
         REAL(KIND=r_dp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_dp), DIMENSION(nx/2+1,ny,nz), INTENT(OUT), TARGET :: in
      END SUBROUTINE            

      SUBROUTINE sfftw_plan_dft_c2r_1d(plan,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: n, flags
         REAL(KIND=r_sp), DIMENSION(n), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_sp), DIMENSION(n/2+1), INTENT(OUT), TARGET :: in
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_c2r_2d(plan,nx,ny,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, flags
         REAL(KIND=r_sp), DIMENSION(nx,ny), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_sp), DIMENSION(nx/2+1,ny), INTENT(OUT), TARGET :: in
      END SUBROUTINE
      SUBROUTINE sfftw_plan_dft_c2r_3d(plan,nx,ny,nz,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: nx, ny, nz, flags
         REAL(KIND=r_sp), DIMENSION(nx,ny,nz), INTENT(IN), TARGET :: out
         COMPLEX(KIND=r_sp), DIMENSION(nx/2+1,ny,nz), INTENT(OUT), TARGET :: in
      END SUBROUTINE            
   END INTERFACE
   
   INTERFACE FFTW_PlanR2C
      SUBROUTINE dfftw_plan_dft_r2c(plan,rank,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         REAL(KIND=r_dp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: in
         COMPLEX(KIND=r_dp), DIMENSION((n(1)/2+1)*PRODUCT(n(2:))), INTENT(IN), TARGET :: out
      END SUBROUTINE   
      SUBROUTINE sfftw_plan_dft_r2c(plan,rank,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         REAL(KIND=r_sp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: in
         COMPLEX(KIND=r_sp), DIMENSION((n(1)/2+1)*PRODUCT(n(2:))), INTENT(IN), TARGET :: out
      END SUBROUTINE   
   END INTERFACE
      
   INTERFACE FFTW_PlanC2R   
      SUBROUTINE dfftw_plan_dft_c2r(plan,rank,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         REAL(KIND=r_dp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: out
         COMPLEX(KIND=r_dp), DIMENSION((n(1)/2+1)*PRODUCT(n(2:))), INTENT(IN), TARGET :: in
      END SUBROUTINE   
      SUBROUTINE sfftw_plan_dft_c2r(plan,rank,n,in,out,flags)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
         INTEGER, INTENT(IN) :: rank, flags
         INTEGER, DIMENSION(rank), INTENT(IN) :: n
         REAL(KIND=r_sp), DIMENSION(PRODUCT(n)), INTENT(OUT), TARGET :: out
         COMPLEX(KIND=r_sp), DIMENSION((n(1)/2+1)*PRODUCT(n(2:))), INTENT(IN), TARGET :: in
      END SUBROUTINE   
   END INTERFACE
   
   INTERFACE FFTW_Execute
      SUBROUTINE dfftw_execute(plan)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
      END SUBROUTINE
   END INTERFACE
   INTERFACE FFTW_DestroyPlan  
      SUBROUTINE dfftw_destroy_plan(plan)
      USE Precision
         INTEGER(KIND=i_dp), INTENT(INOUT) :: plan
      END SUBROUTINE
   END INTERFACE
   
   INTEGER FFTW_R2HC
   PARAMETER (FFTW_R2HC=0)
   INTEGER FFTW_HC2R
   PARAMETER (FFTW_HC2R=1)
   INTEGER FFTW_DHT
   PARAMETER (FFTW_DHT=2)
   INTEGER FFTW_REDFT00
   PARAMETER (FFTW_REDFT00=3)
   INTEGER FFTW_REDFT01
   PARAMETER (FFTW_REDFT01=4)
   INTEGER FFTW_REDFT10
   PARAMETER (FFTW_REDFT10=5)
   INTEGER FFTW_REDFT11
   PARAMETER (FFTW_REDFT11=6)
   INTEGER FFTW_RODFT00
   PARAMETER (FFTW_RODFT00=7)
   INTEGER FFTW_RODFT01
   PARAMETER (FFTW_RODFT01=8)
   INTEGER FFTW_RODFT10
   PARAMETER (FFTW_RODFT10=9)
   INTEGER FFTW_RODFT11
   PARAMETER (FFTW_RODFT11=10)
   INTEGER FFTW_FORWARD
   PARAMETER (FFTW_FORWARD=-1)
   INTEGER FFTW_BACKWARD
   PARAMETER (FFTW_BACKWARD=+1)
   INTEGER FFTW_MEASURE
   PARAMETER (FFTW_MEASURE=0)
   INTEGER FFTW_DESTROY_INPUT
   PARAMETER (FFTW_DESTROY_INPUT=1)
   INTEGER FFTW_UNALIGNED
   PARAMETER (FFTW_UNALIGNED=2)
   INTEGER FFTW_CONSERVE_MEMORY
   PARAMETER (FFTW_CONSERVE_MEMORY=4)
   INTEGER FFTW_EXHAUSTIVE
   PARAMETER (FFTW_EXHAUSTIVE=8)
   INTEGER FFTW_PRESERVE_INPUT
   PARAMETER (FFTW_PRESERVE_INPUT=16)
   INTEGER FFTW_PATIENT
   PARAMETER (FFTW_PATIENT=32)
   INTEGER FFTW_ESTIMATE
   PARAMETER (FFTW_ESTIMATE=64)
   INTEGER FFTW_ESTIMATE_PATIENT
   PARAMETER (FFTW_ESTIMATE_PATIENT=128)
   INTEGER FFTW_BELIEVE_PCOST
   PARAMETER (FFTW_BELIEVE_PCOST=256)
   INTEGER FFTW_DFT_R2HC_ICKY
   PARAMETER (FFTW_DFT_R2HC_ICKY=512)
   INTEGER FFTW_NONTHREADED_ICKY
   PARAMETER (FFTW_NONTHREADED_ICKY=1024)
   INTEGER FFTW_NO_BUFFERING
   PARAMETER (FFTW_NO_BUFFERING=2048)
   INTEGER FFTW_NO_INDIRECT_OP
   PARAMETER (FFTW_NO_INDIRECT_OP=4096)
   INTEGER FFTW_ALLOW_LARGE_GENERIC
   PARAMETER (FFTW_ALLOW_LARGE_GENERIC=8192)
   INTEGER FFTW_NO_RANK_SPLITS
   PARAMETER (FFTW_NO_RANK_SPLITS=16384)
   INTEGER FFTW_NO_VRANK_SPLITS
   PARAMETER (FFTW_NO_VRANK_SPLITS=32768)
   INTEGER FFTW_NO_VRECURSE
   PARAMETER (FFTW_NO_VRECURSE=65536)
   INTEGER FFTW_NO_SIMD
   PARAMETER (FFTW_NO_SIMD=131072)
   
END MODULE FFTW
! EOF
