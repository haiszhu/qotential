#ifndef BIESOLVER_LFORTRAN_ARRAY_H
#define BIESOLVER_LFORTRAN_ARRAY_H

#include <stdbool.h>
#include <complex.h>
#include <stddef.h>
#include <stdint.h>

struct dimension_descriptor {
    int32_t lower_bound;
    int32_t length;
    int32_t stride;
};

struct r64 {
    double *data;
    struct dimension_descriptor dims[32];
    int32_t n_dims;
    int32_t offset;
    bool is_allocated;
};

struct c64 {
    double _Complex *data;
    struct dimension_descriptor dims[32];
    int32_t n_dims;
    int32_t offset;
    bool is_allocated;
};

struct i64 {
    int64_t *data;
    struct dimension_descriptor dims[32];
    int32_t n_dims;
    int32_t offset;
    bool is_allocated;
};

_Static_assert(offsetof(struct r64, dims) == sizeof(double *),
               "unexpected LFortran descriptor dims offset");
_Static_assert(offsetof(struct r64, n_dims) ==
                   sizeof(double *) + 32 * sizeof(struct dimension_descriptor),
               "unexpected LFortran descriptor rank offset");
_Static_assert(offsetof(struct r64, offset) == offsetof(struct r64, n_dims) + 4,
               "unexpected LFortran descriptor offset field");
_Static_assert(offsetof(struct r64, is_allocated) == offsetof(struct r64, offset) + 4,
               "unexpected LFortran descriptor allocation flag");
#if UINTPTR_MAX == UINT64_MAX
_Static_assert(sizeof(struct r64) == 408, "unexpected 64-bit LFortran descriptor size");
#elif UINTPTR_MAX == UINT32_MAX
_Static_assert(sizeof(struct r64) == 400, "unexpected wasm32 LFortran descriptor size");
#else
#error unsupported pointer size
#endif

_Static_assert(sizeof(struct c64) == sizeof(struct r64),
               "unexpected complex descriptor size");
_Static_assert(sizeof(struct i64) == sizeof(struct r64),
               "unexpected integer descriptor size");

int biesolver_rm_to_cm_f64(const double *src, double *dst,
                           size_t rows, size_t cols);
int biesolver_cm_to_rm_f64(const double *src, double *dst,
                           size_t rows, size_t cols);

int64_t biesolver_charts_rowmajor(
    int64_t mp, int64_t np, int64_t p, struct r64 *x0, struct r64 *d,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, double restol, int64_t cap, struct r64 *cl);

int biesolver_geo_charts_rowmajor(
    struct r64 *cl, int64_t nchart, int64_t mp, int64_t np, int64_t p,
    struct r64 *x0, struct r64 *d, struct r64 *uvs, struct r64 *wts,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, struct r64 *sx, struct r64 *snx, struct r64 *sw,
    struct r64 *rts, struct r64 *rps);

int biesolver_uv2x_charts_rowmajor(
    struct r64 *cl, int64_t nchart, int64_t mp, int64_t np, int64_t p,
    struct r64 *x0, int64_t itri, int64_t nuv, struct r64 *uv,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, struct r64 *x);

#endif
