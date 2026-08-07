# qotential/test/Makefile
#
# Throughput benchmark for rrq_r64, ported from
# rrq-legacy/test/stellarator/stellarator.f.  Same driver, same data, the
# only change being rrq -> rrq_r64, so the pps figures are comparable.
#
#   make            build build/stellarator_bench
#   make run        build and run case 1
#   make run N=2    build and run cases 1..2
#   make run N=2 ISIMD=128   same, on the C SIMD kernel (0 / 128 / 256 / 512)
#   make clean
#
# The .dat files are read from rrq-legacy and are not copied here: they run
# from 13 MB (order 4) to 888 MB (order 14), about 1.6 GB in total.
#
#   case  order  mp   np       N        .dat
#     1     4    24   72     44280      13 MB
#     2     6    24   72     93030      28 MB
#     3     8    36  108    358776     110 MB
#     4    10    36  108    548130     168 MB
#     5    12    48  144   1381536     423 MB
#     6    14    60  180               888 MB
#
# The package libraries must be built first:
#   ( cd ../external/LineQuaaadrature && make lib )
#   ( cd ../external/QuatApproximation && make lib )
#   ( cd ..                             && make lib )

ROOT   := $(abspath $(dir $(abspath $(lastword $(MAKEFILE_LIST))))/..)
QA_DIR ?= $(ROOT)/external/QuatApproximation
LQ_DIR ?= $(ROOT)/external/LineQuaaadrature
BLD    := $(ROOT)/build

# Same compiler pick as the package Makefiles -- see the note there.
ifneq ($(filter default environment,$(origin FC)),)
  FC := $(shell for c in gfortran-16 gfortran-15 gfortran-14 gfortran-13 gfortran; do \
                  command -v $$c >/dev/null 2>&1 && { echo $$c; exit 0; }; done; echo gfortran)
endif

FFLAGS := -g -O3 -fPIC -ffp-contract=off -fno-unsafe-math-optimizations \
          -fallow-argument-mismatch -std=legacy -w -fdefault-integer-8 \
          -frecursive -march=native -funroll-loops -fopenmp \
          -I$(BLD) -I$(QA_DIR)/build -I$(LQ_DIR)/build

UNAME_S := $(shell uname -s)
ARCH    := $(shell uname -m)

# Homebrew prefix: /opt/homebrew on Apple Silicon, /usr/local on Intel.
BREW := /opt/homebrew
ifeq ($(UNAME_S),Darwin)
  ifneq ($(ARCH), arm64)
    BREW := /usr/local
  endif
endif

ifeq ($(UNAME_S),Darwin)
  LAPACK_LIBS := -L$(BREW)/opt/openblas-singlethread/lib -lopenblas
  HDF5_LIBS   := -L$(BREW)/opt/hdf5/lib -lhdf5
else
  LAPACK_LIBS := -lopenblas
  HDF5_LIBS   := -lhdf5
endif

LIBS := $(BLD)/libqotential.a \
        $(QA_DIR)/build/libQuatApproximation.a \
        $(LQ_DIR)/build/libLineQuaaadrature.a \
        $(HDF5_LIBS) $(LAPACK_LIBS) -lgfortran -lquadmath -lm

BIN   := $(BLD)/stellarator_bench
SGRF_SRC := stellarator_grf.f90
W7X_SRC := ../web/stellarator/fortran/w7x_modes_mod.f90
W7X_OBJ := $(BLD)/w7x_modes_mod.o
SGRF_CORE_SRC := ../web/stellarator/fortran/stellarator_grf_core_mod.f90
SGRF_CORE_OBJ := $(BLD)/stellarator_grf_core_mod.o
SGRF_BIN := $(BLD)/stellarator_grf
DIRECT_SRC := ../web/stellarator/native/stellarator_direct_main.f90
DIRECT_BIN := $(BLD)/stellarator_direct
GEO_C    := ../utils/stellarator_geo_mex.c
FMM3D_DIR ?= $(HOME)/git/FMM3D
FMM3D_LIB := $(FMM3D_DIR)/lib-static/libfmm3d.a
SGRF_BLAS ?= -L$(BREW)/opt/openblas-threaded/lib -lopenblas
N     ?= 1
ISIMD ?= 0
DATDIR ?= ./

.PHONY: all run sgrf sgrf_direct clean

all: $(BIN)

$(BIN): stellarator_bench.f90 $(BLD)/libqotential.a
	$(FC) $(FFLAGS) $< $(LIBS) -o $@
	@echo "built $@"

sgrf: $(SGRF_BIN)

$(W7X_OBJ): $(W7X_SRC) $(BLD)/libqotential.a
	$(FC) $(FFLAGS) -J$(BLD) -c $< -o $@

$(SGRF_CORE_OBJ): $(SGRF_CORE_SRC) $(W7X_OBJ) $(BLD)/libqotential.a
	$(FC) $(FFLAGS) -cpp -J$(BLD) -c $< -o $@

$(SGRF_BIN): $(SGRF_CORE_OBJ) $(SGRF_SRC) $(GEO_C) $(BLD)/libqotential.a
	$(CC) -O3 -c $(GEO_C) -o $(BLD)/stellarator_geo.o
	$(FC) $(FFLAGS) $(SGRF_SRC) $(SGRF_CORE_OBJ) $(W7X_OBJ) \
	  $(BLD)/stellarator_geo.o $(FMM3D_LIB) \
	  $(BLD)/libqotential.a $(QA_DIR)/build/libQuatApproximation.a \
	  $(LQ_DIR)/build/libLineQuaaadrature.a \
	  $(HDF5_LIBS) $(SGRF_BLAS) -lgfortran -lquadmath -lm -o $@
	@echo "built $@"

sgrf_direct: $(DIRECT_BIN)

$(DIRECT_BIN): $(SGRF_CORE_OBJ) $(W7X_OBJ) $(DIRECT_SRC) $(GEO_C) $(BLD)/libqotential.a
	$(CC) -O3 -c $(GEO_C) -o $(BLD)/stellarator_geo.o
	$(FC) $(FFLAGS) $(DIRECT_SRC) $(SGRF_CORE_OBJ) $(W7X_OBJ) \
	  $(BLD)/stellarator_geo.o $(FMM3D_LIB) $(BLD)/libqotential.a \
	  $(QA_DIR)/build/libQuatApproximation.a $(LQ_DIR)/build/libLineQuaaadrature.a \
	  $(HDF5_LIBS) $(SGRF_BLAS) -lgfortran -lquadmath -lm -o $@
	@echo "built $@"

run: $(BIN)
	cd $(ROOT) && OMP_NUM_THREADS=1 ./build/stellarator_bench $(N) $(ISIMD) $(DATDIR)

clean:
	rm -f $(BIN) $(SGRF_BIN) $(DIRECT_BIN) $(SGRF_CORE_OBJ) $(W7X_OBJ) \
	  $(BLD)/stellarator_grf_core_mod.mod $(BLD)/w7x_modes_mod.mod
