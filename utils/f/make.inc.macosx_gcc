# System-dependent settings for fortran lib

# Change these if you want to use different compilers:
CC = gcc-10
CXX = g++-10
FC = gfortran-10
CFLAGS = -march=native
FFLAGS = -fPIC -w -fallow-argument-mismatch -O3 -march=native -funroll-loops
#CXXFLAGS := $(CFLAGS) -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3

# Change this to location of MATLAB's command line mex executable:
# (For a Mac system, this is something of the form
#  /Applications/MATLAB_R2014a.app/bin/mex )
# For Mac if you want to bypass mex as the compiler, see */makefile.nomex
#FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
#MFLAGS += -L${FDIR} -I/usr/local/Cellar/eigen/3.4.0_1/include/eigen3
MEX = $(shell ls -d /Applications/MATLAB_R20**.app)/bin/mex
# Instead for Octave use: (don't forget to make clean when switch)
#MEX = mkoctave --mex

# Change this to the location of your MWrap executable,
# either an absolute path, or a path relative to
# each of the four demo subdirectories:
MWRAP = ~/mwrap/mwrap

OMPFLAGS = -fopenmp
OMPLIBS = -L/usr/local/lib -lgomp

# FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
# MFLAGS += -L${FDIR}
MFLAGS = -L/usr/local/lib/gcc/10