# makefile for interfacing to a Fortran library
# 
include make.inc

# This is the task run if "make" is done with no arguments ("mexfile" is not
# a file, just a name for a compilation task given below):
default: mexfile

# How to compile your library function: (-fPIC is needed for linking to MEX)
moments.o: moments.f makefile
	$(FC) -fPIC $(FFLAGS) --openmp -c moments.f
harmonics.o: harmonics.f makefile
	$(FC) -fPIC $(FFLAGS) -c harmonics.f
kernels.o: kernels.f makefile
	$(FC) -fPIC $(FFLAGS) --openmp -c kernels.f

# Tell mwrap to make the interface and the name of the final MEX file:
# Note: -c99complex assumes usual Fortran complex numbers (see MWrap docs)
specialquad.c: specialquad.mw makefile
	$(MWRAP) -c99complex -mex specialquad -mb -list specialquad.mw
	$(MWRAP) -c99complex -mex specialquad -c specialquad.c specialquad.mw

# How to make MATLAB/MEX compile the interface (whatever libs needed):
mexfile: specialquad.c moments.o harmonics.o kernels.o makefile
	$(MEX) specialquad.c moments.o harmonics.o kernels.o -largeArrayDims $(MFLAGS) -lgfortran -lm -lgomp

# Move file to bin
mvfile:
	mv -f *.mex* ../bin/
	mv -f *_mex.m ../bin/

# Remove the MEX interface, MATLAB caller, and MEX file:
clean:
	rm -f moments_mex.m specialquad.c *.o *.mex* *~ fort.*
