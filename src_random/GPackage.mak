myf90sources += NURNGs.f90
myf90sources += ParallelRNG.f90
myf90sources += Precision.f90
myf90sources += Random.f90
myf90sources += RNG.f90

# Note that the RNGs are already included in the Parallel version of BoxLib, but not fParallel
# They are extracted here in C and made separate so there is no inconsistencies
mycsources += RNGs.c

# For fParallel-based codes:
f90sources += $(myf90sources) 
csources += $(mycsources)


