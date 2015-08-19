# BOXLIB_HOME defines the directory in which we will find all the BoxLib code
# If you set BOXLIB_HOME as an environment variable, this line will be ignored
BOXLIB_HOME := ../../BoxLib

NDEBUG    := t
MPI       := t
OMP       :=
PROF      :=
COMP      := gfortran
CCOMP     := gcc
MKVERBOSE := t

include $(BOXLIB_HOME)/Tools/F_mk/GMakedefs.mak

VPATH_LOCATIONS += ..

include ../src_multiSpec/GPackage.mak
VPATH_LOCATIONS += ../src_multiSpec

include ../src_gmres/GPackage.mak
VPATH_LOCATIONS += ../src_gmres

include ../src_common/GPackage.mak
VPATH_LOCATIONS += ../src_common

include $(BOXLIB_HOME)/Src/F_BaseLib/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/F_BaseLib

include $(BOXLIB_HOME)/Src/LinearSolvers/F_MG/GPackage.mak
VPATH_LOCATIONS += $(BOXLIB_HOME)/Src/LinearSolvers/F_MG

LAPACK95_DIR = ../src_multiSpec/LAPACK95

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
