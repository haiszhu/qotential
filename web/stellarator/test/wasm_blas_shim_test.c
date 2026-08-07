#include "../native/lfortran_array.h"

#include <complex.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef __EMSCRIPTEN__
#include <pthread.h>
#endif

void biesolver_dgemm(char *ta, char *tb, int64_t m, int64_t n, int64_t k,
                     double alpha, struct r64 *a, int64_t lda,
                     struct r64 *b, int64_t ldb, double beta,
                     struct r64 *c, int64_t ldc);
void biesolver_zgemm(char *ta, char *tb, int64_t m, int64_t n, int64_t k,
                     double _Complex alpha, struct c64 *a, int64_t lda,
                     struct c64 *b, int64_t ldb, double _Complex beta,
                     struct c64 *c, int64_t ldc);
void biesolver_dgesv(int64_t n, int64_t nrhs, struct r64 *a, int64_t lda,
                     struct i64 *ipiv, struct r64 *b, int64_t ldb,
                     int64_t *info);

static void require(int condition, const char *message)
{
    if (!condition) {
        fprintf(stderr, "WASM BLAS shim test failed: %s\n", message);
        exit(1);
    }
}

static struct r64 r64_matrix(double *data, int32_t rows, int32_t cols)
{
    struct r64 value;
    memset(&value, 0, sizeof(value));
    value.data = data;
    value.n_dims = 2;
    value.dims[0] = (struct dimension_descriptor){1, rows, cols};
    value.dims[1] = (struct dimension_descriptor){1, cols, 1};
    value.is_allocated = true;
    return value;
}

static struct c64 c64_matrix(double _Complex *data, int32_t rows, int32_t cols)
{
    struct c64 value;
    memset(&value, 0, sizeof(value));
    value.data = data;
    value.n_dims = 2;
    value.dims[0] = (struct dimension_descriptor){1, rows, cols};
    value.dims[1] = (struct dimension_descriptor){1, cols, 1};
    value.is_allocated = true;
    return value;
}

static struct i64 i64_vector(int64_t *data, int32_t length)
{
    struct i64 value;
    memset(&value, 0, sizeof(value));
    value.data = data;
    value.n_dims = 1;
    value.dims[0] = (struct dimension_descriptor){1, length, 1};
    value.is_allocated = true;
    return value;
}

static int close_real(double actual, double expected)
{
    double scale = fmax(1.0, fabs(expected));
    return fabs(actual - expected) <= 2.0e-14 * scale;
}

static int close_complex(double _Complex actual, double _Complex expected)
{
    double scale = fmax(1.0, cabs(expected));
    return cabs(actual - expected) <= 2.0e-14 * scale;
}

static void test_dgemm_nn(void)
{
    double a_data[6] = {1,2,3, 4,5,6};
    double b_data[12] = {7,8,9,10, 11,12,13,14, 15,16,17,18};
    double c_data[8] = {1,2,3,4, 5,6,7,8};
    const double expected[8] = {92,99,106,113, 213.75,232,250.25,268.5};
    struct r64 a = r64_matrix(a_data, 2, 3);
    struct r64 b = r64_matrix(b_data, 3, 4);
    struct r64 c = r64_matrix(c_data, 2, 4);
    size_t i;
    biesolver_dgemm("N", "N", 2, 4, 3, 1.25, &a, 2, &b, 3,
                     -0.5, &c, 2);
    for (i = 0; i < 8; ++i)
        require(close_real(c_data[i], expected[i]), "non-square dgemm NN mismatch");
}

static void test_dgemm_padded_descriptors(void)
{
    double a_backing[10] = {1,2,3,91,92, 4,5,6,93,94};
    double b_backing[18] = {
        7,8,9,10,81,82, 11,12,13,14,83,84, 15,16,17,18,85,86
    };
    double c_backing[14] = {1,2,3,4,71,72,73, 5,6,7,8,74,75,76};
    const double expected[8] = {92,99,106,113, 213.75,232,250.25,268.5};
    struct r64 a = r64_matrix(a_backing, 2, 5);
    struct r64 b = r64_matrix(b_backing, 3, 6);
    struct r64 c = r64_matrix(c_backing, 2, 7);
    size_t i, j;
    biesolver_dgemm("N", "N", 2, 4, 3, 1.25, &a, 2, &b, 3,
                     -0.5, &c, 2);
    for (i = 0; i < 2; ++i)
        for (j = 0; j < 4; ++j)
            require(close_real(c_backing[i * 7 + j], expected[i * 4 + j]),
                    "padded dgemm descriptor mismatch");
    require(c_backing[4] == 71 && c_backing[5] == 72 && c_backing[6] == 73 &&
            c_backing[11] == 74 && c_backing[12] == 75 && c_backing[13] == 76,
            "padded dgemm overwrote columns outside logical C");
}

static void test_zgemm_tn(void)
{
    double _Complex a_data[6] = {
        1+I, 2-I, -1+2*I,
        3, 0.5+I, 2-2*I
    };
    double _Complex b_data[6] = {1,2+I,-1, 0.5-I,3,2*I};
    double _Complex c_backing[36];
    const double _Complex expected[9] = {
        3-1.5*I, 11+4*I, 0.5+6.5*I,
        5.25+I, 9+5.5*I, -1+5*I,
        1.5+2.5*I, 6+I, 9.5+6.5*I
    };
    struct c64 a = c64_matrix(a_data, 2, 3);
    struct c64 b = c64_matrix(b_data, 2, 3);
    struct c64 c;
    size_t i, j;
    for (i = 0; i < 36; ++i) c_backing[i] = -99-99*I;
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            c_backing[i * 12 + j * 4] = (double)(i * 3 + j + 1) * (1+I);
    c = c64_matrix(c_backing, 3, 3);
    c.dims[0].stride = 12;
    c.dims[1].stride = 4;
    biesolver_zgemm("T", "N", 3, 3, 2, 1.0, &a, 2, &b, 2,
                     0.5, &c, 3);
    for (i = 0; i < 3; ++i)
        for (j = 0; j < 3; ++j)
            require(close_complex(c_backing[i * 12 + j * 4], expected[i * 3 + j]),
                    "strided zgemm TN mismatch");
    require(c_backing[1] == -99-99*I && c_backing[35] == -99-99*I,
            "strided zgemm overwrote an adjacent slice");
}

static void test_dgesv_two_rhs(void)
{
    double a_data[9] = {3,1,-1, 2,4,1, -1,2,5};
    double b_data[6] = {-1,8.5, 1,4, 12,-11};
    const double expected[6] = {1,2, -1,0.5, 3,-2};
    int64_t pivots[3] = {0}, info = -999;
    struct r64 a = r64_matrix(a_data, 3, 3);
    struct r64 b = r64_matrix(b_data, 3, 2);
    struct i64 ipiv = i64_vector(pivots, 3);
    size_t i;
    biesolver_dgesv(3, 2, &a, 3, &ipiv, &b, 3, &info);
    require(info == 0, "dgesv reported failure");
    for (i = 0; i < 6; ++i)
        require(close_real(b_data[i], expected[i]), "dgesv two-RHS mismatch");
    for (i = 0; i < 3; ++i)
        require(pivots[i] >= 1 && pivots[i] <= 3, "dgesv pivot is not 1-based");
}

static void run_numeric_tests(void)
{
    test_dgemm_nn();
    test_dgemm_padded_descriptors();
    test_zgemm_tn();
    test_dgesv_two_rhs();
}

#ifndef __EMSCRIPTEN__
static void *thread_main(void *unused)
{
    (void)unused;
    run_numeric_tests();
    return NULL;
}
#endif

int main(void)
{
    int repeat;
    for (repeat = 0; repeat < 2; ++repeat) run_numeric_tests();
#ifndef __EMSCRIPTEN__
    {
        pthread_t threads[4];
        int i;
        for (i = 0; i < 4; ++i)
            require(pthread_create(&threads[i], NULL, thread_main, NULL) == 0,
                    "pthread_create failed");
        for (i = 0; i < 4; ++i)
            require(pthread_join(threads[i], NULL) == 0, "pthread_join failed");
    }
#endif
    puts("WASM_BLAS_SHIM_OK");
    return 0;
}
