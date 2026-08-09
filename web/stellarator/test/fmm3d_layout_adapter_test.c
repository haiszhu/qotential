#include "../native/fmm3d_layout_adapter.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int allocation_count;
static int fail_allocation_at;
static int fmm_call_count;

void *biesolver_test_malloc(size_t size)
{
    ++allocation_count;
    if (allocation_count == fail_allocation_at) return NULL;
    return malloc(size);
}

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "FMM3D layout adapter test failed: %s\n", message);
        exit(1);
    }
}

static struct r64 matrix(double *data, int32_t rows, int32_t cols)
{
    struct r64 result;
    memset(&result, 0, sizeof(result));
    result.data = data;
    result.dims[0] = (struct dimension_descriptor){1, rows, cols};
    result.dims[1] = (struct dimension_descriptor){1, cols, 1};
    result.n_dims = 2;
    result.is_allocated = true;
    return result;
}

static struct r64 vector(double *data, int32_t length)
{
    struct r64 result;
    memset(&result, 0, sizeof(result));
    result.data = data;
    result.dims[0] = (struct dimension_descriptor){1, length, 1};
    result.n_dims = 1;
    result.is_allocated = true;
    return result;
}

static void interleave(const double *planar, double *column_major, int64_t n)
{
    int64_t i, j;
    for (j = 0; j < n; ++j)
        for (i = 0; i < 3; ++i)
            column_major[3*j+i] = planar[i*n+j];
}

static void direct_oracle(int64_t ns, const double *source,
                          const double *charge, const double *dipvec,
                          int64_t nt, const double *target, double *potential)
{
    const long double inv4pi =
        7.957747154594766788444188168626e-2L;
    int64_t i, j;
    for (i = 0; i < nt; ++i) {
        long double value = 0.0L;
        for (j = 0; j < ns; ++j) {
            long double dx = (long double)target[i] - source[j];
            long double dy = (long double)target[nt+i] - source[ns+j];
            long double dz = (long double)target[2*nt+i] - source[2*ns+j];
            long double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 == 0.0L) continue;
            long double rinv = 1.0L/sqrtl(r2);
            long double dot = dx*dipvec[j] + dy*dipvec[ns+j] +
                              dz*dipvec[2*ns+j];
            value += inv4pi*(charge[j]*rinv + dot*rinv/r2);
        }
        potential[i] = (double)value;
    }
}

static double scaled_error(const double *actual, const double *expected,
                           int64_t n)
{
    double numerator = 0.0, denominator = 1.0;
    int64_t i;
    for (i = 0; i < n; ++i) {
        numerator = fmax(numerator, fabs(actual[i]-expected[i]));
        denominator = fmax(denominator, fabs(expected[i]));
    }
    return numerator/denominator;
}

#ifndef TEST_ALLOCATION_FAILURE_ONLY
void lfmm3d_t_cd_p_(double *, int64_t *, double *, double *, double *,
                     int64_t *, double *, double *, int64_t *);

static void test_layout_and_distinct_targets(void)
{
    enum { NS = 5, NT = 4 };
    double source_data[3*NS] = {
        -0.7, 0.2, 1.1, -1.3, 0.45,
         0.4,-0.9, 0.3,  1.2,-0.15,
         0.8, 0.6,-1.0,  0.1, 1.35
    };
    double target_data[3*NT] = {
         0.1,-0.8, 0.75, 1.4,
        -0.2, 0.55,1.1, -0.6,
         1.6,-1.2, 0.05,0.7
    };
    double charge_data[NS] = {0.7,-1.2,0.35,1.1,-0.4};
    double dipole_data[3*NS] = {
         0.2,-0.1, 0.45, 0.3,-0.5,
        -0.4, 0.6,-0.25, 0.8, 0.15,
         0.9, 0.05,-0.7,-0.2, 0.4
    };
    double actual[NT], expected[NT], wrong[NT];
    double source_cm[3*NS], target_cm[3*NT], dipole_cm[3*NS];
    struct r64 source = matrix(source_data, 3, NS);
    struct r64 target = matrix(target_data, 3, NT);
    struct r64 charge = vector(charge_data, NS);
    struct r64 dipole = matrix(dipole_data, 3, NS);
    struct r64 output = vector(actual, NT);
    int64_t ier = -1, status = -1, raw_ier = -1;
    double elapsed = -1.0, eps = 1.0e-12;

    direct_oracle(NS, source_data, charge_data, dipole_data,
                  NT, target_data, expected);
    biesolver_lfmm3d_t_cd_p_rowmajor(
        eps, NS, &source, &charge, &dipole, NT, &target, &output,
        &ier, &status, &elapsed);
    require(status == 0, "adapter rejected valid nonsymmetric descriptors");
    require(ier == 0, "FMM3D reported an error");
    require(elapsed >= 0.0, "elapsed time is negative");
    require(scaled_error(actual, expected, NT) <= 2.0e-11,
            "FMM differs from the long-double direct oracle");

    interleave(source_data, source_cm, NS);
    interleave(target_data, target_cm, NT);
    interleave(dipole_data, dipole_cm, NS);
    memset(wrong, 0, sizeof(wrong));
    /* Deliberately pass planar storage to the column-major ABI. */
    lfmm3d_t_cd_p_(&eps, (int64_t[]){NS}, source_data, charge_data,
                   dipole_data, (int64_t[]){NT}, target_data, wrong, &raw_ier);
    require(raw_ier == 0, "negative-control FMM call failed");
    require(scaled_error(wrong, expected, NT) > 1.0e-4,
            "wrong-layout negative control did not detect the ABI seam");

    /* The converted arrays are also checked explicitly against a raw call. */
    memset(wrong, 0, sizeof(wrong));
    lfmm3d_t_cd_p_(&eps, (int64_t[]){NS}, source_cm, charge_data,
                   dipole_cm, (int64_t[]){NT}, target_cm, wrong, &raw_ier);
    require(raw_ier == 0 && scaled_error(wrong, actual, NT) <= 2.0e-14,
            "adapter does not match the raw column-major ABI");
    puts("FMM3D_LAYOUT_PARITY_OK");
    puts("FMM3D_LAYOUT_NEGATIVE_CONTROL_OK");
}

static void test_same_shape_distinct_target(void)
{
    enum { N = 5 };
    double source_data[3*N] = {
        -1.0,-0.5,0.1,0.8,1.4, 0.2,1.1,-0.7,0.4,-1.2,
         0.9,-0.3,1.5,-0.8,0.25
    };
    double target_data[3*N] = {
         0.3,1.7,-1.1,0.55,-0.2, -1.4,0.2,0.85,-0.6,1.3,
        -0.9,0.65,0.4,1.8,-1.5
    };
    double charge_data[N] = {1.0,-0.2,0.4,-0.7,0.9};
    double dipole_data[3*N] = {
        0.1,0.3,-0.5,0.7,-0.9, 0.2,-0.4,0.6,-0.8,1.0,
        0.9,-0.7,0.5,-0.3,0.1
    };
    double actual[N], expected[N];
    struct r64 source = matrix(source_data, 3, N);
    struct r64 target = matrix(target_data, 3, N);
    struct r64 charge = vector(charge_data, N);
    struct r64 dipole = matrix(dipole_data, 3, N);
    struct r64 output = vector(actual, N);
    int64_t ier = -1, status = -1;
    double elapsed = -1.0;

    direct_oracle(N, source_data, charge_data, dipole_data,
                  N, target_data, expected);
    biesolver_lfmm3d_t_cd_p_rowmajor(
        1.0e-12, N, &source, &charge, &dipole, N, &target, &output,
        &ier, &status, &elapsed);
    require(status == 0 && ier == 0, "same-shape distinct-target call failed");
    require(scaled_error(actual, expected, N) <= 2.0e-11,
            "same-shape target was incorrectly aliased to source");
}

static void test_tree_path_source_target_alias(void)
{
    enum { N = 2048 };
    const double pi = 3.141592653589793238462643383279502884;
    const double golden = 0.6180339887498948482;
    double *source_data = calloc(3*N, sizeof(double));
    double *charge_data = calloc(N, sizeof(double));
    double *dipole_data = calloc(3*N, sizeof(double));
    double *actual = calloc(N, sizeof(double));
    double *expected = calloc(N, sizeof(double));
    struct r64 source, charge, dipole, output;
    int64_t ier = -1, status = -1;
    double elapsed = -1.0;
    int i;
    require(source_data && charge_data && dipole_data && actual && expected,
            "tree-path fixture allocation failed");
    for (i = 0; i < N; ++i) {
        double z = -0.95 + 1.9*((double)i + 0.5)/N;
        double theta = 2*pi*fmod(golden*i, 1.0);
        double radius = sqrt(1-z*z);
        source_data[i] = radius*cos(theta);
        source_data[N+i] = radius*sin(theta);
        source_data[2*N+i] = z;
        charge_data[i] = sin(0.017*i) + 0.25*cos(0.031*i);
        dipole_data[i] = 0.1*cos(0.013*i);
        dipole_data[N+i] = -0.07*sin(0.019*i);
        dipole_data[2*N+i] = 0.05*cos(0.023*i);
    }
    source = matrix(source_data, 3, N);
    charge = vector(charge_data, N);
    dipole = matrix(dipole_data, 3, N);
    output = vector(actual, N);
    direct_oracle(N, source_data, charge_data, dipole_data,
                  N, source_data, expected);
    biesolver_lfmm3d_t_cd_p_rowmajor(
        1.0e-9, N, &source, &charge, &dipole, N, &source, &output,
        &ier, &status, &elapsed);
    require(status == 0 && ier == 0, "tree-path alias call failed");
    for (i = 0; i < N; ++i)
        require(isfinite(actual[i]), "tree-path alias produced non-finite output");
    require(scaled_error(actual, expected, N) <= 2.0e-8,
            "tree-path alias differs from direct oracle");
    free(expected); free(actual); free(dipole_data); free(charge_data);
    free(source_data);
}

int main(void)
{
    test_layout_and_distinct_targets();
    test_same_shape_distinct_target();
    test_tree_path_source_target_alias();
    return 0;
}
#else
void biesolver_test_fmm_call(double *eps, int64_t *ns, double *source,
                             double *charge, double *dipole, int64_t *nt,
                             double *target, double *output, int64_t *ier)
{
    (void)eps; (void)ns; (void)source; (void)charge; (void)dipole;
    (void)nt; (void)target; (void)output;
    ++fmm_call_count;
    *ier = 0;
}

int main(void)
{
    enum { N = 5 };
    double source_data[3*N] = {0}, target_data[3*N] = {0};
    double charge_data[N] = {0}, dipole_data[3*N] = {0}, output_data[N];
    struct r64 source = matrix(source_data, 3, N);
    struct r64 target = matrix(target_data, 3, N);
    struct r64 charge = vector(charge_data, N);
    struct r64 dipole = matrix(dipole_data, 3, N);
    struct r64 output = vector(output_data, N);
    int64_t ier = -1, status = -1;
    double elapsed = -1.0;
    int i;
    for (i = 0; i < N; ++i) output_data[i] = 99.0;
    allocation_count = 0;
    fail_allocation_at = 2;
    fmm_call_count = 0;
    biesolver_lfmm3d_t_cd_p_rowmajor(
        1.0e-6, N, &source, &charge, &dipole, N, &target, &output,
        &ier, &status, &elapsed);
    require(status == 1, "allocation failure did not set adapter status");
    require(ier == 0 && elapsed == 0.0, "allocation failure leaked FMM status");
    require(fmm_call_count == 0, "FMM was called after adapter allocation failure");
    for (i = 0; i < N; ++i)
        require(output_data[i] == 0.0, "allocation failure did not zero output");
    puts("FMM3D_LAYOUT_ALLOCATION_FAILURE_OK");
    return 0;
}
#endif
