include $(HYDROLIB_HOME)/GPackageCommon.mak

# We do not include ParallelRNG.f90 here since it requires all of fParallel/BoxLib

libraries += -L/usr/lib64/ -lfftw3 -lfftw3f -lm # Dynamic FFTW
#libraries += ${HPC_LIBS_64}/static/libfftw3*.a -lm # Static (manual)
#libraries += -Wl,-Bstatic -L/usr/lib64/ -lfftw3 -lfftw3f -Wl,-Bdynamic -lm # Static FFTW (using ld)
