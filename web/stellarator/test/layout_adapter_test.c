#include "../native/lfortran_array.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

static size_t rm(size_t i, size_t j, size_t n2) { return i * n2 + j; }
static size_t cm(size_t i, size_t j, size_t n1) { return i + n1 * j; }

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "layout adapter test failed: %s\n", message);
        exit(1);
    }
}

static struct r64 descriptor(double *data, int32_t rank,
                             int32_t n0, int32_t n1)
{
    struct r64 value;
    memset(&value, 0, sizeof(value));
    value.data = data;
    value.n_dims = rank;
    value.dims[0].lower_bound = 1;
    value.dims[0].length = n0;
    value.dims[0].stride = 1;
    if (rank == 2) {
        value.dims[1].lower_bound = 1;
        value.dims[1].length = n1;
        value.dims[1].stride = n0;
    }
    value.is_allocated = true;
    return value;
}

static void test_nonsquare_conversions(void)
{
    const double row_major[6] = {11, 12, 13, 21, 22, 23};
    double column_major[6] = {0};
    double roundtrip[6] = {0};
    size_t i, j;

    require(biesolver_rm_to_cm_f64(row_major, column_major, 2, 3) == 0,
            "2x3 rm->cm failed");
    for (i = 0; i < 2; ++i)
        for (j = 0; j < 3; ++j)
            require(column_major[cm(i, j, 2)] == row_major[rm(i, j, 3)],
                    "2x3 rm->cm value mismatch");
    require(biesolver_cm_to_rm_f64(column_major, roundtrip, 2, 3) == 0,
            "2x3 cm->rm failed");
    require(memcmp(row_major, roundtrip, sizeof(row_major)) == 0,
            "2x3 conversion did not round-trip");
    require(biesolver_rm_to_cm_f64(row_major, (double *)row_major, 2, 3) != 0,
            "in-place rm->cm must be rejected");
    require(biesolver_cm_to_rm_f64(row_major, (double *)row_major, 2, 3) != 0,
            "in-place cm->rm must be rejected");

    {
        double bands[15], interleaved[15], bands_again[15];
        for (i = 0; i < 3; ++i)
            for (j = 0; j < 5; ++j)
                bands[rm(i, j, 5)] = 100.0 * (double)(i + 1) + (double)j;
        require(biesolver_rm_to_cm_f64(bands, interleaved, 3, 5) == 0,
                "3x5 vector rm->cm failed");
        require(biesolver_cm_to_rm_f64(interleaved, bands_again, 3, 5) == 0,
                "3x5 vector cm->rm failed");
        require(memcmp(bands, bands_again, sizeof(bands)) == 0,
                "3x5 coordinate bands did not round-trip");
    }
}

static void gauss4_differentiation(const double *nodes, double *d_rm)
{
    double weights[4];
    size_t i, j, k;
    for (i = 0; i < 4; ++i) {
        weights[i] = 1.0;
        for (j = 0; j < 4; ++j)
            if (j != i) weights[i] /= nodes[i] - nodes[j];
    }
    for (i = 0; i < 4; ++i) {
        double diagonal = 0.0;
        for (j = 0; j < 4; ++j) {
            if (i == j) continue;
            k = rm(i, j, 4);
            d_rm[k] = weights[j] / (weights[i] * (nodes[i] - nodes[j]));
            diagonal -= d_rm[k];
        }
        d_rm[rm(i, i, 4)] = diagonal;
    }
}

static void test_production_boundaries(void)
{
    const i8 mp = 4, np = 12, p = 4, nfp = 0, nmode = 0, cap = 4096;
    const i8 hdim = 10;
    const double restol = 1.0e300;
    double x0_data[4] = {
        -0.8611363115940525752, -0.3399810435848562648,
         0.3399810435848562648,  0.8611363115940525752
    };
    double d_rm[16] = {0}, d_cm[16] = {0};
    double *cl_wrapped = calloc((size_t)(6 * cap), sizeof(double));
    double *cl_direct = calloc((size_t)(6 * cap), sizeof(double));
    double dummy = 0.0;
    struct r64 x0 = descriptor(x0_data, 1, 4, 0);
    struct r64 d = descriptor(d_rm, 2, 4, 4);
    struct r64 mn = descriptor(&dummy, 1, 1, 0);
    struct r64 rc = descriptor(&dummy, 1, 1, 0);
    struct r64 zs = descriptor(&dummy, 1, 1, 0);
    struct r64 cl = descriptor(cl_wrapped, 1, (int32_t)(6 * cap), 0);
    i8 nchart_wrapped, nchart_direct, nsrc;
    size_t vector_count;

    require(cl_wrapped != NULL && cl_direct != NULL, "chart allocation failed");
    gauss4_differentiation(x0_data, d_rm);
    require(biesolver_rm_to_cm_f64(d_rm, d_cm, 4, 4) == 0,
            "D conversion failed");
    nchart_wrapped = biesolver_charts_rowmajor(
        mp, np, p, &x0, &d, nfp, nmode, &mn, &rc, &zs, restol, cap, &cl);
    nchart_direct = stellarator_charts(
        &mp, &np, &p, x0_data, d_cm, &nfp, &nmode, &dummy, &dummy, &dummy,
        &restol, &cap, cl_direct);
    require(nchart_wrapped > 0 && nchart_wrapped == nchart_direct,
            "chart wrapper count differs from production");
    require(memcmp(cl_wrapped, cl_direct,
                   (size_t)(6 * nchart_direct) * sizeof(double)) == 0,
            "chart wrapper output differs from production");

    nsrc = 2 * nchart_direct * hdim;
    vector_count = (size_t)(3 * nsrc);
    {
        double uvs_rm[20] = {
            0.0, 0.25, 0.5, 0.75, 0.0, 0.25, 0.5, 0.0, 0.25, 0.0,
            0.0, 0.0,  0.0, 0.0,  0.25,0.25, 0.25,0.5, 0.5,  0.75
        };
        double uvs_cm[20], wts[10];
        double *sx_rm = calloc(vector_count, sizeof(double));
        double *snx_rm = calloc(vector_count, sizeof(double));
        double *rts_rm = calloc(vector_count, sizeof(double));
        double *rps_rm = calloc(vector_count, sizeof(double));
        double *sw_rm = calloc((size_t)nsrc, sizeof(double));
        double *sx_cm = calloc(vector_count, sizeof(double));
        double *snx_cm = calloc(vector_count, sizeof(double));
        double *rts_cm = calloc(vector_count, sizeof(double));
        double *rps_cm = calloc(vector_count, sizeof(double));
        double *sw_direct = calloc((size_t)nsrc, sizeof(double));
        double *expected = calloc(vector_count, sizeof(double));
        struct r64 uvs = descriptor(uvs_rm, 2, 2, 10);
        struct r64 wts_desc = descriptor(wts, 1, 10, 0);
        struct r64 sx = descriptor(sx_rm, 2, 3, (int32_t)nsrc);
        struct r64 snx = descriptor(snx_rm, 2, 3, (int32_t)nsrc);
        struct r64 sw = descriptor(sw_rm, 1, (int32_t)nsrc, 0);
        struct r64 rts = descriptor(rts_rm, 2, 3, (int32_t)nsrc);
        struct r64 rps = descriptor(rps_rm, 2, 3, (int32_t)nsrc);
        size_t i;

        require(sx_rm && snx_rm && rts_rm && rps_rm && sw_rm && sx_cm &&
                snx_cm && rts_cm && rps_cm && sw_direct && expected,
                "geometry allocation failed");
        for (i = 0; i < 10; ++i) wts[i] = 0.05 + 0.001 * (double)i;
        biesolver_rm_to_cm_f64(uvs_rm, uvs_cm, 2, 10);
        require(biesolver_geo_charts_rowmajor(
                    &cl, nchart_wrapped, mp, np, p, &x0, &d, &uvs, &wts_desc,
                    nfp, nmode, &mn, &rc, &zs, &sx, &snx, &sw, &rts, &rps) == 0,
                "geometry wrapper failed");
        stellarator_geo_charts(
            cl_direct, &nchart_direct, &mp, &np, &p, x0_data, d_cm,
            uvs_cm, wts, &nfp, &nmode, &dummy, &dummy, &dummy,
            sx_cm, snx_cm, sw_direct, rts_cm, rps_cm);
#define REQUIRE_VECTOR_MATCH(cm_values, rm_values, label) do { \
        biesolver_cm_to_rm_f64((cm_values), expected, 3, (size_t)nsrc); \
        require(memcmp(expected, (rm_values), vector_count * sizeof(double)) == 0, \
                (label)); \
    } while (0)
        REQUIRE_VECTOR_MATCH(sx_cm, sx_rm, "sx layout differs");
        REQUIRE_VECTOR_MATCH(snx_cm, snx_rm, "snx layout differs");
        REQUIRE_VECTOR_MATCH(rts_cm, rts_rm, "rts layout differs");
        REQUIRE_VECTOR_MATCH(rps_cm, rps_rm, "rps layout differs");
#undef REQUIRE_VECTOR_MATCH
        require(memcmp(sw_direct, sw_rm, (size_t)nsrc * sizeof(double)) == 0,
                "sw differs from production");
        free(sx_rm); free(snx_rm); free(rts_rm); free(rps_rm); free(sw_rm);
        free(sx_cm); free(snx_cm); free(rts_cm); free(rps_cm); free(sw_direct);
        free(expected);
    }

    {
        const i8 itri = 1, nuv = 5;
        double uv_rm[10] = {0,0.2,0.4,0.1,0.3, 0,0.1,0.2,0.5,0.3};
        double uv_cm[10], x_rm[15] = {0}, x_cm[15] = {0}, expected[15];
        struct r64 uv = descriptor(uv_rm, 2, 2, 5);
        struct r64 x = descriptor(x_rm, 2, 3, 5);
        biesolver_rm_to_cm_f64(uv_rm, uv_cm, 2, 5);
        require(biesolver_uv2x_charts_rowmajor(
                    &cl, nchart_wrapped, mp, np, p, &x0, itri, nuv, &uv,
                    nfp, nmode, &mn, &rc, &zs, &x) == 0,
                "uv2x wrapper failed");
        stellarator_uv2x_charts(
            cl_direct, &nchart_direct, &mp, &np, &p, x0_data, &itri, &nuv,
            uv_cm, &nfp, &nmode, &dummy, &dummy, &dummy, x_cm);
        biesolver_cm_to_rm_f64(x_cm, expected, 3, 5);
        require(memcmp(expected, x_rm, sizeof(expected)) == 0,
                "uv2x wrapper output differs from production");
    }
    free(cl_wrapped);
    free(cl_direct);
}

int main(void)
{
    test_nonsquare_conversions();
    test_production_boundaries();
    puts("LAYOUT_ADAPTER_OK");
    return 0;
}
