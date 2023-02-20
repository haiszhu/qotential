# makefile overrides
# OS:       macOS
# Compiler: gfortran 9.X
# OpenMP:   enabled
#

CC=/opt/homebrew/opt/llvm/bin/clang
CXX=/opt/homebrew/opt/llvm/bin/clang++
FC=gfortran-12

ifeq ($(PREFIX),)
    FMM_INSTALL_DIR=/Users/hzhu/lib
endif

CFLAGS = -O3 -arch arm64
FFLAGS = $(CFLAGS) -std=legacy
CFLAGS += -I src 

# OpenMP with gcc on OSX needs the following
OMPFLAGS = -fopenmp
OMPLIBS = -lgomp

# MATLAB interface:
FDIR=$$(dirname `gfortran --print-file-name libgfortran.dylib`)
MFLAGS += -L${FDIR}
MEX = /Applications/MATLAB_R20***.app/bin/mex
#LIBS = -lm -lstdc++.6
#MEXLIBS= -lm -lstdc++.6 -lgfortran -ldl

CFLAGS += -I/opt/homebrew/opt/libomp/include
OMPFLAGS = -Xpreprocessor -fopenmp
OMPLIBS = $(shell ls -d /Applications/MATLAB_R2022b_Beta.app)/toolbox/eml/externalDependency/omp/maca64/lib/libomp.dylib


