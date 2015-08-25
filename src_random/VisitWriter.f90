MODULE VisIt_Writer
IMPLICIT NONE
PUBLIC

INTERFACE

!void write_point_mesh(const char *filename, int useBinary, int npts, 
!                      float *pts, int nvars, int *vardim, 
!                      const char * const *varnames, float **vars);
SUBROUTINE WritePointVTKMesh (filename, ub, npts, pts, nvars, vardim, varnames, vars) &
   BIND(C,NAME="write_point_mesh")
   USE ISO_C_BINDING
   TYPE(C_PTR), VALUE :: filename
   !CHARACTER (C_CHAR), DIMENSION (*) :: filename
   INTEGER (C_INT), VALUE :: ub, npts, nvars
   REAL (C_DOUBLE), DIMENSION (3*npts) :: pts
   INTEGER (C_INT), DIMENSION (nvars) :: vardim
   TYPE(C_PTR), DIMENSION (nvars) :: varnames, vars
END SUBROUTINE
   
!void write_regular_mesh(const char *filename, int useBinary, int *dims, 
!                        int nvars, int *vardim, int *centering,
!                        const char * const *varnames, float **vars);
SUBROUTINE WriteRegularVTKMesh (filename, ub, dims, nvars, vardim, &
   centering, varnames, vars) BIND(C,NAME="write_regular_mesh")
   USE ISO_C_BINDING
   TYPE(C_PTR), VALUE :: filename
   !CHARACTER (C_CHAR), DIMENSION (*) :: filename
   INTEGER (C_INT), VALUE :: ub, nvars
   INTEGER (C_INT), DIMENSION (3) :: dims
   INTEGER (C_INT), DIMENSION (nvars) :: vardim, centering
   TYPE(C_PTR), DIMENSION (nvars) :: varnames, vars
END SUBROUTINE

!void write_rectilinear_mesh(const char *filename, int useBinary, 
!                            int *dims, float *x, float *y, float *z, 
!                            int nvars, int *vardim, int *centering, 
!                            const char * const *varnames, float **vars);
SUBROUTINE WriteRectilinearVTKMesh (filename, ub, dims, x, y, z, nvars, &
   vardim, centering, varnames, vars) BIND(C,NAME="write_rectilinear_mesh")
   USE ISO_C_BINDING
   TYPE(C_PTR), VALUE :: filename
   !CHARACTER (C_CHAR), DIMENSION (*) :: filename
   INTEGER (C_INT), VALUE :: ub, nvars
   INTEGER (C_INT), DIMENSION (3) :: dims
   REAL (C_DOUBLE), DIMENSION (dims(1)) :: x
   REAL (C_DOUBLE), DIMENSION (dims(2)) :: y
   REAL (C_DOUBLE), DIMENSION (dims(3)) :: z
   INTEGER (C_INT), DIMENSION (nvars) :: vardim, centering
   TYPE(C_PTR), DIMENSION (nvars) :: varnames, vars
END SUBROUTINE

END INTERFACE

END MODULE
