#include "lfortran_array.h"

#include <stdint.h>
#include <stdlib.h>

typedef int64_t i8;

int64_t stellarator_charts(const i8 *mp, const i8 *np, const i8 *p,
                           const double *x0, const double *d,
                           const i8 *nfp, const i8 *nmode, const double *mn,
                           const double *rc, const double *zs,
                           const double *restol, const i8 *cap, double *cl);
void stellarator_geo_charts(const double *cl, const i8 *nchart,
                            const i8 *mp, const i8 *np, const i8 *p,
                            const double *x0, const double *d,
                            const double *uvs, const double *wts,
                            const i8 *nfp, const i8 *nmode, const double *mn,
                            const double *rc, const double *zs,
                            double *sx, double *snx, double *sw,
                            double *rts, double *rps);
void stellarator_uv2x_charts(const double *cl, const i8 *nchart,
                             const i8 *mp, const i8 *np, const i8 *p,
                             const double *x0, const i8 *itri, const i8 *nuv,
                             const double *uv, const i8 *nfp,
                             const i8 *nmode, const double *mn,
                             const double *rc, const double *zs, double *x);
i8 stellarator_geo_ntri(const i8 *mp, const i8 *np, const i8 *p,
                        const double *x0);
void stellarator_geo(const i8 *mp, const i8 *np, const i8 *p,
                     const double *x0, const double *d,
                     const double *uvs, const double *wts,
                     double *sx, double *snx, double *sw,
                     double *rts, double *rps);
void stellarator_geo_uv2x(const i8 *mp, const i8 *np, const i8 *p,
                          const double *x0, const i8 *itri, const i8 *nuv,
                          const double *uv, double *x);

static int checked_product(size_t a, size_t b, size_t *result)
{
    if (a != 0 && b > SIZE_MAX / a) return 0;
    *result = a * b;
    return 1;
}

static int valid_extent(int64_t value)
{
    return value >= 0 && (uint64_t)value <= (uint64_t)SIZE_MAX;
}

static int checked_i64_product(int64_t a, int64_t b, int64_t *result)
{
    if (a < 0 || b < 0 || (a != 0 && b > INT64_MAX / a)) return 0;
    *result = a * b;
    return 1;
}

static double *descriptor_data(struct r64 *array)
{
    if (array == NULL || array->data == NULL || array->offset < 0) return NULL;
    return array->data + (size_t)array->offset;
}

static int descriptor_has_shape(const struct r64 *array, int rank,
                                int64_t n0, int64_t n1)
{
    if (array == NULL || array->data == NULL || array->offset < 0 ||
        array->n_dims != rank || !valid_extent(n0) ||
        (rank == 2 && !valid_extent(n1))) return 0;
    if ((int64_t)array->dims[0].length < n0) return 0;
    if (rank == 2 && (int64_t)array->dims[1].length < n1) return 0;
    return 1;
}

int biesolver_rm_to_cm_f64(const double *src, double *dst,
                           size_t rows, size_t cols)
{
    size_t count, i, j;
    if (!checked_product(rows, cols, &count)) return 2;
    if (count == 0) return 0;
    if (src == NULL || dst == NULL || src == dst) return 1;
    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j)
            dst[i + rows * j] = src[i * cols + j];
    return 0;
}

int biesolver_cm_to_rm_f64(const double *src, double *dst,
                           size_t rows, size_t cols)
{
    size_t count, i, j;
    if (!checked_product(rows, cols, &count)) return 2;
    if (count == 0) return 0;
    if (src == NULL || dst == NULL || src == dst) return 1;
    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j)
            dst[i * cols + j] = src[i + rows * j];
    return 0;
}

int64_t biesolver_charts_rowmajor(
    int64_t mp, int64_t np, int64_t p, struct r64 *x0, struct r64 *d,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, double restol, int64_t cap, struct r64 *cl)
{
    size_t d_count;
    int64_t cl_count, mn_count;
    double *d_cm;
    int64_t result;

    if (mp <= 0 || np <= 0 || p <= 0 || nmode < 0 || cap <= 0 ||
        !checked_i64_product(6, cap, &cl_count) ||
        !checked_i64_product(2, nmode, &mn_count) ||
        !descriptor_has_shape(x0, 1, p, 0) ||
        !descriptor_has_shape(d, 2, p, p) ||
        !descriptor_has_shape(cl, 1, cl_count, 0) ||
        (nmode > 0 && (!descriptor_has_shape(mn, 1, mn_count, 0) ||
                       !descriptor_has_shape(rc, 1, nmode, 0) ||
                       !descriptor_has_shape(zs, 1, nmode, 0)))) return -2;
    if (!valid_extent(p) || !checked_product((size_t)p, (size_t)p, &d_count) ||
        d_count > SIZE_MAX / sizeof(double)) return -3;
    d_cm = (double *)malloc(d_count * sizeof(double));
    if (d_cm == NULL) return -4;
    if (biesolver_rm_to_cm_f64(descriptor_data(d), d_cm,
                               (size_t)p, (size_t)p) != 0) {
        free(d_cm);
        return -2;
    }
    result = stellarator_charts(&mp, &np, &p, descriptor_data(x0), d_cm,
                                &nfp, &nmode, descriptor_data(mn),
                                descriptor_data(rc), descriptor_data(zs),
                                &restol, &cap, descriptor_data(cl));
    free(d_cm);
    return result;
}

int biesolver_geo_charts_rowmajor(
    struct r64 *cl, int64_t nchart, int64_t mp, int64_t np, int64_t p,
    struct r64 *x0, struct r64 *d, struct r64 *uvs, struct r64 *wts,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, struct r64 *sx, struct r64 *snx, struct r64 *sw,
    struct r64 *rts, struct r64 *rps)
{
    size_t hdim, nsrc, d_count, uv_count, vector_count, total_count;
    int64_t p_plus_one, hdim_i64, ntri_i64, nsrc_i64, cl_count, mn_count;
    double *workspace, *d_cm, *uvs_cm, *sx_cm, *snx_cm, *rts_cm, *rps_cm;

    if (nchart <= 0 || mp <= 0 || np <= 0 || p <= 0 || nmode < 0 ||
        p == INT64_MAX || !valid_extent(nchart) || !valid_extent(p)) return 1;
    p_plus_one = p + 1;
    if (!checked_i64_product(p, p_plus_one, &hdim_i64)) return 2;
    hdim_i64 /= 2;
    if (!checked_i64_product(2, nchart, &ntri_i64) ||
        !checked_i64_product(ntri_i64, hdim_i64, &nsrc_i64) ||
        !checked_i64_product(6, nchart, &cl_count) ||
        !checked_i64_product(2, nmode, &mn_count) ||
        !valid_extent(hdim_i64) || !valid_extent(nsrc_i64)) return 2;
    hdim = (size_t)hdim_i64;
    nsrc = (size_t)nsrc_i64;
    if (
        !checked_product((size_t)p, (size_t)p, &d_count) ||
        !checked_product(2, hdim, &uv_count) ||
        !checked_product(3, nsrc, &vector_count) ||
        !checked_product(4, vector_count, &total_count) ||
        total_count > SIZE_MAX - d_count ||
        total_count + d_count > SIZE_MAX - uv_count ||
        total_count + d_count + uv_count > SIZE_MAX / sizeof(double)) return 2;
    if (!descriptor_has_shape(cl, 1, cl_count, 0) ||
        !descriptor_has_shape(x0, 1, p, 0) ||
        !descriptor_has_shape(d, 2, p, p) ||
        !descriptor_has_shape(uvs, 2, 2, hdim_i64) ||
        !descriptor_has_shape(wts, 1, hdim_i64, 0) ||
        !descriptor_has_shape(sx, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(snx, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(sw, 1, nsrc_i64, 0) ||
        !descriptor_has_shape(rts, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(rps, 2, 3, nsrc_i64) ||
        (nmode > 0 && (!descriptor_has_shape(mn, 1, mn_count, 0) ||
                       !descriptor_has_shape(rc, 1, nmode, 0) ||
                       !descriptor_has_shape(zs, 1, nmode, 0)))) return 1;

    workspace = (double *)malloc((d_count + uv_count + total_count) * sizeof(double));
    if (workspace == NULL) return 3;
    d_cm = workspace;
    uvs_cm = d_cm + d_count;
    sx_cm = uvs_cm + uv_count;
    snx_cm = sx_cm + vector_count;
    rts_cm = snx_cm + vector_count;
    rps_cm = rts_cm + vector_count;
    biesolver_rm_to_cm_f64(descriptor_data(d), d_cm, (size_t)p, (size_t)p);
    biesolver_rm_to_cm_f64(descriptor_data(uvs), uvs_cm, 2, hdim);
    stellarator_geo_charts(
        descriptor_data(cl), &nchart, &mp, &np, &p, descriptor_data(x0),
        d_cm, uvs_cm, descriptor_data(wts), &nfp, &nmode,
        descriptor_data(mn), descriptor_data(rc), descriptor_data(zs),
        sx_cm, snx_cm, descriptor_data(sw), rts_cm, rps_cm);
    biesolver_cm_to_rm_f64(sx_cm, descriptor_data(sx), 3, nsrc);
    biesolver_cm_to_rm_f64(snx_cm, descriptor_data(snx), 3, nsrc);
    biesolver_cm_to_rm_f64(rts_cm, descriptor_data(rts), 3, nsrc);
    biesolver_cm_to_rm_f64(rps_cm, descriptor_data(rps), 3, nsrc);
    free(workspace);
    return 0;
}

int biesolver_uv2x_charts_rowmajor(
    struct r64 *cl, int64_t nchart, int64_t mp, int64_t np, int64_t p,
    struct r64 *x0, int64_t itri, int64_t nuv, struct r64 *uv,
    int64_t nfp, int64_t nmode, struct r64 *mn, struct r64 *rc,
    struct r64 *zs, struct r64 *x)
{
    size_t uv_count, x_count, total_count;
    int64_t ntri_i64, cl_count, mn_count;
    double *workspace, *uv_cm, *x_cm;

    if (nchart <= 0 || mp <= 0 || np <= 0 || p <= 0 || itri <= 0 ||
        nuv <= 0 || nmode < 0 ||
        !checked_i64_product(2, nchart, &ntri_i64) || itri > ntri_i64 ||
        !checked_i64_product(6, nchart, &cl_count) ||
        !checked_i64_product(2, nmode, &mn_count) || !valid_extent(nuv) ||
        !checked_product(2, (size_t)nuv, &uv_count) ||
        !checked_product(3, (size_t)nuv, &x_count) ||
        uv_count > SIZE_MAX - x_count ||
        (total_count = uv_count + x_count) > SIZE_MAX / sizeof(double)) return 1;
    if (!descriptor_has_shape(cl, 1, cl_count, 0) ||
        !descriptor_has_shape(x0, 1, p, 0) ||
        !descriptor_has_shape(uv, 2, 2, nuv) ||
        !descriptor_has_shape(x, 2, 3, nuv) ||
        (nmode > 0 && (!descriptor_has_shape(mn, 1, mn_count, 0) ||
                       !descriptor_has_shape(rc, 1, nmode, 0) ||
                       !descriptor_has_shape(zs, 1, nmode, 0)))) return 1;
    workspace = (double *)malloc(total_count * sizeof(double));
    if (workspace == NULL) return 3;
    uv_cm = workspace;
    x_cm = uv_cm + uv_count;
    biesolver_rm_to_cm_f64(descriptor_data(uv), uv_cm, 2, (size_t)nuv);
    stellarator_uv2x_charts(
        descriptor_data(cl), &nchart, &mp, &np, &p, descriptor_data(x0),
        &itri, &nuv, uv_cm, &nfp, &nmode, descriptor_data(mn),
        descriptor_data(rc), descriptor_data(zs), x_cm);
    biesolver_cm_to_rm_f64(x_cm, descriptor_data(x), 3, (size_t)nuv);
    free(workspace);
    return 0;
}

/* Built-in-surface entry points (stellarator_geo_ntri / _geo / _geo_uv2x),
   same seam as the charts adapters above: the generated C hands over
   row-major descriptors, the production C computes in column-major. */

int64_t biesolver_geo_ntri_rowmajor(int64_t mp, int64_t np, int64_t p,
                                    struct r64 *x0)
{
    if (mp <= 0 || np <= 0 || p <= 0 ||
        !descriptor_has_shape(x0, 1, p, 0)) return -2;
    return stellarator_geo_ntri(&mp, &np, &p, descriptor_data(x0));
}

int biesolver_geo_rowmajor(
    int64_t ntri, int64_t mp, int64_t np, int64_t p,
    struct r64 *x0, struct r64 *d, struct r64 *uvs, struct r64 *wts,
    struct r64 *sx, struct r64 *snx, struct r64 *sw,
    struct r64 *rts, struct r64 *rps)
{
    size_t hdim, nsrc, d_count, uv_count, vector_count, total_count;
    int64_t p_plus_one, hdim_i64, nsrc_i64;
    double *workspace, *d_cm, *uvs_cm, *sx_cm, *snx_cm, *rts_cm, *rps_cm;

    if (ntri <= 0 || mp <= 0 || np <= 0 || p <= 0 ||
        p == INT64_MAX || !valid_extent(ntri) || !valid_extent(p)) return 1;
    p_plus_one = p + 1;
    if (!checked_i64_product(p, p_plus_one, &hdim_i64)) return 2;
    hdim_i64 /= 2;
    if (!checked_i64_product(ntri, hdim_i64, &nsrc_i64) ||
        !valid_extent(hdim_i64) || !valid_extent(nsrc_i64)) return 2;
    hdim = (size_t)hdim_i64;
    nsrc = (size_t)nsrc_i64;
    if (!checked_product((size_t)p, (size_t)p, &d_count) ||
        !checked_product(2, hdim, &uv_count) ||
        !checked_product(3, nsrc, &vector_count) ||
        !checked_product(4, vector_count, &total_count) ||
        total_count > SIZE_MAX - d_count ||
        total_count + d_count > SIZE_MAX - uv_count ||
        total_count + d_count + uv_count > SIZE_MAX / sizeof(double)) return 2;
    if (!descriptor_has_shape(x0, 1, p, 0) ||
        !descriptor_has_shape(d, 2, p, p) ||
        !descriptor_has_shape(uvs, 2, 2, hdim_i64) ||
        !descriptor_has_shape(wts, 1, hdim_i64, 0) ||
        !descriptor_has_shape(sx, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(snx, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(sw, 1, nsrc_i64, 0) ||
        !descriptor_has_shape(rts, 2, 3, nsrc_i64) ||
        !descriptor_has_shape(rps, 2, 3, nsrc_i64)) return 1;

    workspace = (double *)malloc((d_count + uv_count + total_count) * sizeof(double));
    if (workspace == NULL) return 3;
    d_cm = workspace;
    uvs_cm = d_cm + d_count;
    sx_cm = uvs_cm + uv_count;
    snx_cm = sx_cm + vector_count;
    rts_cm = snx_cm + vector_count;
    rps_cm = rts_cm + vector_count;
    biesolver_rm_to_cm_f64(descriptor_data(d), d_cm, (size_t)p, (size_t)p);
    biesolver_rm_to_cm_f64(descriptor_data(uvs), uvs_cm, 2, hdim);
    stellarator_geo(&mp, &np, &p, descriptor_data(x0), d_cm, uvs_cm,
                    descriptor_data(wts), sx_cm, snx_cm, descriptor_data(sw),
                    rts_cm, rps_cm);
    biesolver_cm_to_rm_f64(sx_cm, descriptor_data(sx), 3, nsrc);
    biesolver_cm_to_rm_f64(snx_cm, descriptor_data(snx), 3, nsrc);
    biesolver_cm_to_rm_f64(rts_cm, descriptor_data(rts), 3, nsrc);
    biesolver_cm_to_rm_f64(rps_cm, descriptor_data(rps), 3, nsrc);
    free(workspace);
    return 0;
}

int biesolver_geo_uv2x_rowmajor(
    int64_t mp, int64_t np, int64_t p, struct r64 *x0,
    int64_t itri, int64_t nuv, struct r64 *uv, struct r64 *x)
{
    size_t uv_count, x_count, total_count;
    double *workspace, *uv_cm, *x_cm;

    if (mp <= 0 || np <= 0 || p <= 0 || itri <= 0 || nuv <= 0 ||
        !valid_extent(nuv) ||
        !checked_product(2, (size_t)nuv, &uv_count) ||
        !checked_product(3, (size_t)nuv, &x_count) ||
        uv_count > SIZE_MAX - x_count ||
        (total_count = uv_count + x_count) > SIZE_MAX / sizeof(double)) return 1;
    if (!descriptor_has_shape(x0, 1, p, 0) ||
        !descriptor_has_shape(uv, 2, 2, nuv) ||
        !descriptor_has_shape(x, 2, 3, nuv)) return 1;
    workspace = (double *)malloc(total_count * sizeof(double));
    if (workspace == NULL) return 3;
    uv_cm = workspace;
    x_cm = uv_cm + uv_count;
    biesolver_rm_to_cm_f64(descriptor_data(uv), uv_cm, 2, (size_t)nuv);
    stellarator_geo_uv2x(&mp, &np, &p, descriptor_data(x0), &itri, &nuv,
                         uv_cm, x_cm);
    biesolver_cm_to_rm_f64(x_cm, descriptor_data(x), 3, (size_t)nuv);
    free(workspace);
    return 0;
}
