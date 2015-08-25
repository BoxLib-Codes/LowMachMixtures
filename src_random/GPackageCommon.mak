vpath %.f90 $(HYDROLIB_HOME)
vpath %.c $(HYDROLIB_HOME)
vpath %.h $(HYDROLIB_HOME)
INCLUDE_LOCATIONS += $(HYDROLIB_HOME)

myf90sources += Precision.f90
myf90sources += Random.f90
myf90sources += NURNGs.f90
myf90sources += RNG.f90
myf90sources += FFTW.f90
myf90sources += VisitWriter.f90
myf90sources += HydroGridModule.f90
myf90sources += HydroGridCInterface.f90

# Note that the RNGs are already included in the Parallel version of BoxLib, but not fParallel
# They are extracted here in C and made separate so there is no inconsistencies
mycsources += RNGs.c visit_writer.c
mycheaders += RNGs.h visit_writer.h

# For fParallel-based codes:
f90sources += $(myf90sources) 
csources += $(mycsources)
cheaders += $(mycheaders) # This does not actually seem to exist in fParallel?!?

# For Parallel-based codes:
f90EXE_sources += $(myf90sources)
cEXE_sources   += $(mycsources)
cEXE_headers   += $(mycheaders)
