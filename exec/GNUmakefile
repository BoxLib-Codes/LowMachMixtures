# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME ?= ../../BoxLib
LOWMACH_HOME ?= ../

NDEBUG    := t
MPI       := t
OMP       :=
PROF      :=
MKVERBOSE := t

# Choose the compiler (make sure the LAPACK95 interface is compiled by the same compiler)
COMP      := gfortran
CCOMP     := gcc

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

VPATH_LOCATIONS += $(LOWMACH_HOME)

include $(LOWMACH_HOME)/src_multiSpec/GPackage.mak
VPATH_LOCATIONS += $(LOWMACH_HOME)/src_multiSpec

include $(LOWMACH_HOME)/src_gmres/GPackage.mak
VPATH_LOCATIONS += $(LOWMACH_HOME)/src_gmres

include $(LOWMACH_HOME)/src_common/GPackage.mak
VPATH_LOCATIONS += $(LOWMACH_HOME)/src_common

include $(LOWMACH_HOME)/src_random/GPackage.mak
VPATH_LOCATIONS += $(LOWMACH_HOME)/src_random

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

# Make sure to pre-compile LAPACK95 and have lapack/blas available:
LAPACK95_DIR = $(LOWMACH_HOME)/src_lapack
F90FLAGS += -I$(LAPACK95_DIR)/lapack95_modules
ifeq ($(findstring hopper, $(HOST)), hopper)
  libraries += $(LAPACK95_DIR)/lapack95.a
else ifeq ($(findstring edison, $(HOST)), edison)
  libraries += $(LAPACK95_DIR)/lapack95.a
else ifeq ($(findstring bint01, $(HOST)), bint01)
  libraries += $(LAPACK95_DIR)/lapack95.a
else
  libraries += $(LAPACK95_DIR)/lapack95.a -llapack -lblas
endif

main.$(suf).exe: $(objects) 
	$(LINK.f90) -o main.$(suf).exe $(objects) $(libraries)

include $(BOXLIB_HOME)/Tools/F_mk/GMakerules.mak
