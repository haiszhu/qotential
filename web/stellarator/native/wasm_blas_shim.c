#include "lfortran_array.h"

#include <complex.h>
#include <math.h>
#include <stddef.h>
#include <stdint.h>

static int transpose_kind(const char *flag)
{
    if (flag == NULL) return -1;
    if (*flag == 'N' || *flag == 'n') return 0;
    if (*flag == 'T' || *flag == 't') return 1;
    if (*flag == 'C' || *flag == 'c') return 2;
    return -1;
}

static int valid_r64_matrix(const struct r64 *a, int64_t rows, int64_t cols)
{
    return a != NULL && a->data != NULL && a->n_dims == 2 && a->offset >= 0 &&
           rows >= 0 && cols >= 0 && a->dims[0].length >= rows &&
           a->dims[1].length >= cols && a->dims[0].stride != 0 &&
           a->dims[1].stride != 0;
}

static int valid_c64_matrix(const struct c64 *a, int64_t rows, int64_t cols)
{
    return a != NULL && a->data != NULL && a->n_dims == 2 && a->offset >= 0 &&
           rows >= 0 && cols >= 0 && a->dims[0].length >= rows &&
           a->dims[1].length >= cols && a->dims[0].stride != 0 &&
           a->dims[1].stride != 0;
}

static ptrdiff_t matrix_index(const struct dimension_descriptor *dims,
                              int32_t offset, int64_t row, int64_t col)
{
    return (ptrdiff_t)offset + (ptrdiff_t)row * dims[0].stride +
           (ptrdiff_t)col * dims[1].stride;
}

static double r64_get(const struct r64 *a, int64_t row, int64_t col)
{
    return a->data[matrix_index(a->dims, a->offset, row, col)];
}

static void r64_set(struct r64 *a, int64_t row, int64_t col, double value)
{
    a->data[matrix_index(a->dims, a->offset, row, col)] = value;
}

static double _Complex c64_get(const struct c64 *a, int64_t row, int64_t col)
{
    return a->data[matrix_index(a->dims, a->offset, row, col)];
}

static void c64_set(struct c64 *a, int64_t row, int64_t col,
                    double _Complex value)
{
    a->data[matrix_index(a->dims, a->offset, row, col)] = value;
}

static int valid_leading_dimensions(int ta, int tb, int64_t m, int64_t n,
                                    int64_t k, int64_t lda, int64_t ldb,
                                    int64_t ldc)
{
    int64_t need_a = ta == 0 ? m : k;
    int64_t need_b = tb == 0 ? k : n;
    if (need_a < 1) need_a = 1;
    if (need_b < 1) need_b = 1;
    return lda >= need_a && ldb >= need_b && ldc >= (m > 1 ? m : 1);
}

void biesolver_dgemm(char *ta_flag, char *tb_flag,
                     int64_t m, int64_t n, int64_t k,
                     double alpha, struct r64 *a, int64_t lda,
                     struct r64 *b, int64_t ldb, double beta,
                     struct r64 *c, int64_t ldc)
{
    int ta = transpose_kind(ta_flag), tb = transpose_kind(tb_flag);
    int64_t a_rows, a_cols, b_rows, b_cols, i, j, q;
    if (ta < 0 || tb < 0 || ta == 2 || tb == 2 || m < 0 || n < 0 || k < 0 ||
        !valid_leading_dimensions(ta, tb, m, n, k, lda, ldb, ldc)) return;
    a_rows = ta == 0 ? m : k;
    a_cols = ta == 0 ? k : m;
    b_rows = tb == 0 ? k : n;
    b_cols = tb == 0 ? n : k;
    if (!valid_r64_matrix(a, a_rows, a_cols) ||
        !valid_r64_matrix(b, b_rows, b_cols) ||
        !valid_r64_matrix(c, m, n)) return;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            double sum = 0.0;
            for (q = 0; q < k; ++q) {
                double av = ta == 0 ? r64_get(a, i, q) : r64_get(a, q, i);
                double bv = tb == 0 ? r64_get(b, q, j) : r64_get(b, j, q);
                sum += av * bv;
            }
            double value = alpha * sum;
            if (beta != 0.0) value += beta * r64_get(c, i, j);
            r64_set(c, i, j, value);
        }
    }
}

void biesolver_zgemm(char *ta_flag, char *tb_flag,
                     int64_t m, int64_t n, int64_t k,
                     double _Complex alpha, struct c64 *a, int64_t lda,
                     struct c64 *b, int64_t ldb, double _Complex beta,
                     struct c64 *c, int64_t ldc)
{
    int ta = transpose_kind(ta_flag), tb = transpose_kind(tb_flag);
    int64_t a_rows, a_cols, b_rows, b_cols, i, j, q;
    if (ta < 0 || tb < 0 || m < 0 || n < 0 || k < 0 ||
        !valid_leading_dimensions(ta, tb, m, n, k, lda, ldb, ldc)) return;
    a_rows = ta == 0 ? m : k;
    a_cols = ta == 0 ? k : m;
    b_rows = tb == 0 ? k : n;
    b_cols = tb == 0 ? n : k;
    if (!valid_c64_matrix(a, a_rows, a_cols) ||
        !valid_c64_matrix(b, b_rows, b_cols) ||
        !valid_c64_matrix(c, m, n)) return;
    for (i = 0; i < m; ++i) {
        for (j = 0; j < n; ++j) {
            double _Complex sum = 0.0;
            for (q = 0; q < k; ++q) {
                double _Complex av = ta == 0 ? c64_get(a, i, q) : c64_get(a, q, i);
                double _Complex bv = tb == 0 ? c64_get(b, q, j) : c64_get(b, j, q);
                if (ta == 2) av = conj(av);
                if (tb == 2) bv = conj(bv);
                sum += av * bv;
            }
            double _Complex value = alpha * sum;
            if (creal(beta) != 0.0 || cimag(beta) != 0.0)
                value += beta * c64_get(c, i, j);
            c64_set(c, i, j, value);
        }
    }
}

void biesolver_dgesv(int64_t n, int64_t nrhs, struct r64 *a, int64_t lda,
                     struct i64 *ipiv, struct r64 *b, int64_t ldb,
                     int64_t *info)
{
    int64_t k, i, j, rhs;
    if (info == NULL) return;
    *info = 0;
    if (n < 0) { *info = -1; return; }
    if (nrhs < 0) { *info = -2; return; }
    if (!valid_r64_matrix(a, n, n)) { *info = -3; return; }
    if (lda < (n > 1 ? n : 1)) { *info = -4; return; }
    if (ipiv == NULL || ipiv->data == NULL || ipiv->n_dims != 1 ||
        ipiv->offset < 0 || ipiv->dims[0].length < n ||
        ipiv->dims[0].stride == 0) { *info = -5; return; }
    if (!valid_r64_matrix(b, n, nrhs)) { *info = -6; return; }
    if (ldb < (n > 1 ? n : 1)) { *info = -7; return; }

    for (k = 0; k < n; ++k) {
        int64_t pivot = k;
        double pivot_abs = fabs(r64_get(a, k, k));
        for (i = k + 1; i < n; ++i) {
            double candidate = fabs(r64_get(a, i, k));
            if (candidate > pivot_abs) {
                pivot_abs = candidate;
                pivot = i;
            }
        }
        ipiv->data[(ptrdiff_t)ipiv->offset +
                   (ptrdiff_t)k * ipiv->dims[0].stride] = pivot + 1;
        if (pivot_abs == 0.0) {
            *info = k + 1;
            return;
        }
        if (pivot != k) {
            for (j = 0; j < n; ++j) {
                double tmp = r64_get(a, k, j);
                r64_set(a, k, j, r64_get(a, pivot, j));
                r64_set(a, pivot, j, tmp);
            }
            for (rhs = 0; rhs < nrhs; ++rhs) {
                double tmp = r64_get(b, k, rhs);
                r64_set(b, k, rhs, r64_get(b, pivot, rhs));
                r64_set(b, pivot, rhs, tmp);
            }
        }
        for (i = k + 1; i < n; ++i) {
            double factor = r64_get(a, i, k) / r64_get(a, k, k);
            r64_set(a, i, k, factor);
            for (j = k + 1; j < n; ++j)
                r64_set(a, i, j, r64_get(a, i, j) - factor * r64_get(a, k, j));
            for (rhs = 0; rhs < nrhs; ++rhs)
                r64_set(b, i, rhs,
                        r64_get(b, i, rhs) - factor * r64_get(b, k, rhs));
        }
    }
    for (rhs = 0; rhs < nrhs; ++rhs) {
        for (i = n; i-- > 0;) {
            double value = r64_get(b, i, rhs);
            for (j = i + 1; j < n; ++j)
                value -= r64_get(a, i, j) * r64_get(b, j, rhs);
            r64_set(b, i, rhs, value / r64_get(a, i, i));
        }
    }
}
