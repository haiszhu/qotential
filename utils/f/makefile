# makefile for interfacing to a Fortran library
# 
include make.inc

UNAME := $(shell uname)
ARCH  := $(shell uname -m)

ifeq ($(UNAME), Darwin)
  MATLAB_ROOT := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  ifeq ($(ARCH), arm64)
    MATLAB_ARCH := maca64
    MEX_EXT := mexmaca64
  else
    MATLAB_ARCH := maci64
    MEX_EXT := mexmaci64
  endif
  MATLAB_INC := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                 $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  MEX_LDFLAGS := -bundle -Wl,-undefined,dynamic_lookup
endif

MEX_CC := gcc-15

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
	$(MEX_CC) $(MEX_LDFLAGS) -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  $(MATLAB_INC) \
	  specialquad.c moments.o harmonics.o kernels.o \
	  $(MATLAB_LIBS) \
	  $(MFLAGS) -lgfortran -lm -lgomp \
	  -o specialquad.$(MEX_EXT)

# Move file to bin
mvfile:
	mv -f *.mex* ../bin/
	mv -f *_mex.m ../bin/

# Remove the MEX interface, MATLAB caller, and MEX file:
clean:
	rm -f moments_mex.m specialquad.c *.o *.mex* *~ fort.*
