libraries += -L/usr/lib64 -L$(HYDROLIB_HOME) -lHydroGrid -L${HPC_LIBS_64}/static -lTAU -lDISLIN -lOpenGLUT -lnfft -L/usr/X11R6/lib64 -ltiff -lGLU -lGL -lXmu -lX11 -lXext -lXxf86vm -lXi
# It is useful to link these statically to make a more portable executable:
libraries += /usr/lib64/libfftw3f.a /usr/lib64/libfftw3.a
