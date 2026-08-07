/*
 * C-interoperable BLAS boundary for the LFortran stellarator build.
 *
 * A Fortran BLAS call with CHARACTER arguments has compiler-specific hidden
 * length arguments.  LFortran's bind(C) correctly removes those arguments,
 * so bind(C) must not point directly at the Fortran dgemm_/zgemm_ symbols.
 * Route GEMM through the standardized CBLAS ABI instead.  DGESV has no
 * CHARACTER arguments and is forwarded to the ILP64 Fortran LAPACK symbol.
 */

#include <ctype.h>
#include <stdint.h>

#include <cblas.h>

/* qotential is built with 64-bit default integers, so every size and leading
   dimension arriving here is int64_t.  blasint follows the linked OpenBLAS: it
   is 64-bit only when that build defines OPENBLAS_USE64BITINT.  Against an
   LP64 build the casts below would truncate silently and return wrong results
   for large problems, so make the mismatch a build error instead. */
_Static_assert(sizeof(blasint) == sizeof(int64_t),
               "link an ILP64 BLAS: qotential passes 64-bit integers");

static enum CBLAS_TRANSPOSE
transpose_from_char(char trans)
{
    switch (toupper((unsigned char)trans)) {
    case 'N': return CblasNoTrans;
    case 'T': return CblasTrans;
    case 'C': return CblasConjTrans;
    default:  return CblasNoTrans;
    }
}

void
biesolver_dgemm(char transa, char transb,
                const int64_t *m, const int64_t *n, const int64_t *k,
                const double *alpha, const double *a, const int64_t *lda,
                const double *b, const int64_t *ldb, const double *beta,
                double *c, const int64_t *ldc)
{
    cblas_dgemm(CblasColMajor,
                transpose_from_char(transa), transpose_from_char(transb),
                (blasint)*m, (blasint)*n, (blasint)*k, *alpha,
                a, (blasint)*lda, b, (blasint)*ldb, *beta, c, (blasint)*ldc);
}

void
biesolver_zgemm(char transa, char transb,
                const int64_t *m, const int64_t *n, const int64_t *k,
                const void *alpha, const void *a, const int64_t *lda,
                const void *b, const int64_t *ldb, const void *beta,
                void *c, const int64_t *ldc)
{
    cblas_zgemm(CblasColMajor,
                transpose_from_char(transa), transpose_from_char(transb),
                (blasint)*m, (blasint)*n, (blasint)*k, alpha,
                a, (blasint)*lda, b, (blasint)*ldb, beta, c, (blasint)*ldc);
}

extern void dgesv_(const int64_t *n, const int64_t *nrhs,
                   double *a, const int64_t *lda, int64_t *ipiv,
                   double *b, const int64_t *ldb, int64_t *info);

void
biesolver_dgesv(const int64_t *n, const int64_t *nrhs,
                double *a, const int64_t *lda, int64_t *ipiv,
                double *b, const int64_t *ldb, int64_t *info)
{
    dgesv_(n, nrhs, a, lda, ipiv, b, ldb, info);
}
