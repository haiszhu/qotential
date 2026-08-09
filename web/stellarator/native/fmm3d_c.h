#ifndef BIESOLVER_FMM3D_C_H
#define BIESOLVER_FMM3D_C_H

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef FMM3D_DROP_IN
#define FNAME(x) x##_
#else
#define FNAME(x) x##_c_
#endif

typedef int32_t fint;
typedef int64_t flong;
typedef double _Complex fcomplex;

#define FA2(i, j, ld1) (((j) - 1) * (ld1) + ((i) - 1))
#define FA3(i, j, k, ld1, ld2) \
    ((((k) - 1) * (ld2) + ((j) - 1)) * (ld1) + ((i) - 1))
#define FA4(i, j, k, l, ld1, ld2, ld3) \
    ((((((l) - 1) * (ld3) + ((k) - 1)) * (ld2)) + ((j) - 1)) * \
      (ld1) + ((i) - 1))

#endif
