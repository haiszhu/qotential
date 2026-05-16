# qotential/Makefile
#
# Builds the qotential BIE-solver glue layer:
#   build/libqotential.a              -- Fortran static library
#   matlab/qotential_mex.<EXT>        -- MATLAB-callable mex bundle
#   build/test_lap3d_close            -- pure-Fortran smoke test
#
# Links against two upstream packages (their `make lib` must already
# have been run): QuatApproximation and LineQuaaadrature.  Default
# locations point at the public-repo checkouts; override on the command
# line if your layout differs:
#
#   make mex QA_DIR=/path/to/QuatApproximation \
#            LQ_DIR=/path/to/LineQuaaadrature
#
# Usage:
#   make            -- show targets
#   make lib        -- build build/libqotential.a
#   make mex        -- build matlab/qotential_mex.<EXT>
#   make test       -- build build/test_lap3d_close
#   make clean      -- remove build artifacts

SHELL := /bin/bash

ROOT       := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
ROOT       := $(patsubst %/,%,$(ROOT))
SRC_DIR    := $(ROOT)/src
MATLAB_DIR := $(ROOT)/matlab
TEST_DIR   := $(ROOT)/test
BLD_DIR    := $(ROOT)/build

# ---- upstream package locations (override on command line) ----
QA_DIR ?= /Users/hzhu/Documents/Github/QuatApproximation
LQ_DIR ?= /Users/hzhu/Documents/Github/LineQuaaadrature

QA_LIB := $(QA_DIR)/build/libQuatApproximation.a
LQ_LIB := $(LQ_DIR)/build/libLineQuaaadrature.a
QA_INC := -I$(QA_DIR)/build
LQ_INC := -I$(LQ_DIR)/build

# ---- compilers ----
FC := gfortran-15
CC := gcc-15
MW := ~/mwrap/mwrap
MWFLAGS := -c99complex -i8 -mex

# ---- platform detection ----
UNAME := $(shell uname)
ARCH  := $(shell uname -m)

ifeq ($(UNAME), Darwin)
  MATLAB_ROOT  := $(shell ls -d /Applications/MATLAB_R*.app 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH  := maca64
  MEX_EXT      := mexmaca64
  OPENBLAS_DIR := /opt/homebrew/opt/openblas-singlethread
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmx.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmex.dylib \
                  $(MATLAB_ROOT)/bin/$(MATLAB_ARCH)/libmat.dylib -lm
  MEX_LDFLAGS  := -bundle -Wl,-undefined,dynamic_lookup
  LAPACK_LIBS  := -framework Accelerate
else
  MATLAB_ROOT  := $(shell ls -d /usr/local/MATLAB/R* 2>/dev/null | sort | tail -n1)
  MATLAB_ARCH  := glnxa64
  MEX_EXT      := mexa64
  MATLAB_INC   := -I$(MATLAB_ROOT)/extern/include
  MATLAB_LIBS  := -L$(MATLAB_ROOT)/bin/$(MATLAB_ARCH) -lmx -lmex -lmat -lm
  MEX_LDFLAGS  := -shared
  LAPACK_LIBS  := -L$(OPENBLAS_DIR)/lib -lopenblas
endif

HDF5_ROOT ?= /opt/homebrew/opt/hdf5
HDF5_LIBS := -L$(HDF5_ROOT)/lib -lhdf5

# ---- Fortran flags (identical to both legacy packages) ----
FFLAGS := -g -O3 -fPIC \
           -ffp-contract=off -fno-unsafe-math-optimizations \
           -fallow-argument-mismatch -std=legacy -w \
           -fdefault-integer-8 -frecursive -cpp \
           -march=native -funroll-loops -fopenmp \
           -J$(BLD_DIR) -I$(BLD_DIR) $(QA_INC) $(LQ_INC)

# ---- sources ----
QOT_SOURCES := $(SRC_DIR)/lap3d_close_mod.f90 \
               $(SRC_DIR)/lap3d_close_mex.f90

QOT_OBJECTS := $(patsubst $(SRC_DIR)/%.f90,$(BLD_DIR)/%.o,$(QOT_SOURCES))
LIB         := $(BLD_DIR)/libqotential.a

# ---- mwrap-generated gateway ----
MW_SRC  := $(MATLAB_DIR)/qotential.mw
MEX_C   := $(MATLAB_DIR)/qotential_mex.c
MEX_OUT := $(MATLAB_DIR)/qotential_mex.$(MEX_EXT)

# ---- test ----
TEST_SRC := $(TEST_DIR)/test_lap3d_close.f90
TEST_BIN := $(BLD_DIR)/test_lap3d_close

# ============================================================
.PHONY: all lib mex test clean check-deps

all:
	@echo "qotential -- BIE-solver glue (r64 mex interface)"
	@echo ""
	@echo "  make lib         build $(notdir $(LIB))"
	@echo "  make mex         build $(notdir $(MEX_OUT))"
	@echo "  make test        build $(notdir $(TEST_BIN))"
	@echo "  make clean       remove build artifacts"
	@echo ""
	@echo "Upstream package paths (override on command line):"
	@echo "  QA_DIR = $(QA_DIR)"
	@echo "  LQ_DIR = $(LQ_DIR)"

check-deps:
	@test -f $(QA_LIB) || { echo "ERROR: $(QA_LIB) not found.  Run 'make lib' in $(QA_DIR) first."; exit 1; }
	@test -f $(LQ_LIB) || { echo "ERROR: $(LQ_LIB) not found.  Run 'make lib' in $(LQ_DIR) first."; exit 1; }

lib: $(LIB)

mex: $(MEX_OUT)

test: $(TEST_BIN)

$(BLD_DIR):
	mkdir -p $(BLD_DIR)

# ---- compile Fortran ----
$(BLD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BLD_DIR) check-deps
	$(FC) $(FFLAGS) -c $< -o $@

# module dependency: mex wrapper needs the module first
$(BLD_DIR)/lap3d_close_mex.o: $(BLD_DIR)/lap3d_close_mod.o

# ---- static library ----
$(LIB): $(QOT_OBJECTS)
	ar rcs $@ $^

# ---- mwrap: two-pass generation ----
$(MEX_C): $(MW_SRC) | $(BLD_DIR)
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) qotential_mex -mb -list qotential.mw
	cd $(MATLAB_DIR) && $(MW) $(MWFLAGS) qotential_mex -c qotential_mex.c qotential.mw
	perl -pi -e 's/_{2,}/_/g' $(MEX_C)

# ---- link MEX bundle ----
$(MEX_OUT): $(LIB) $(MEX_C) | check-deps
	$(CC) $(MEX_LDFLAGS) -fPIC \
	  -DMATLAB_MEX_FILE -DMATLAB_DEFAULT_RELEASE=R2018a -DMX_COMPAT_32=0 \
	  $(MATLAB_INC) \
	  $(MEX_C) \
	  -L$(BLD_DIR) -lqotential \
	  $(QA_LIB) $(LQ_LIB) \
	  $(HDF5_LIBS) $(MATLAB_LIBS) $(LAPACK_LIBS) \
	  -lgfortran -lquadmath -lm \
	  -o $(MEX_OUT)

# ---- pure-Fortran test binary ----
$(TEST_BIN): $(LIB) $(TEST_SRC) | $(BLD_DIR)
	$(FC) $(FFLAGS) $(TEST_SRC) \
	  -L$(BLD_DIR) -lqotential \
	  $(QA_LIB) $(LQ_LIB) \
	  $(HDF5_LIBS) $(LAPACK_LIBS) \
	  -lgfortran -lquadmath -lm \
	  -o $(TEST_BIN)

clean:
	rm -rf $(BLD_DIR) $(MEX_C) $(MEX_OUT) $(MATLAB_DIR)/qol_*.m
