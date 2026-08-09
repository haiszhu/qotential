#include "fmm3d_layout_adapter.h"

#include <emscripten/emscripten.h>
#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#ifndef BIESOLVER_FMM3D_MALLOC
#define BIESOLVER_FMM3D_MALLOC malloc
#else
void *BIESOLVER_FMM3D_MALLOC(size_t size);
#endif

#ifndef BIESOLVER_FMM3D_CALL
#define BIESOLVER_FMM3D_CALL lfmm3d_t_cd_p_
#endif
void BIESOLVER_FMM3D_CALL(double *eps, int64_t *nsource, double *source,
                          double *charge, double *dipvec, int64_t *ntarg,
                          double *targ, double *pottarg, int64_t *ier);

static int checked_count(int64_t n, size_t factor, size_t *result)
{
    if (n <= 0 || (uint64_t)n > (uint64_t)(SIZE_MAX/factor)) return 0;
    *result = factor*(size_t)n;
    return *result <= SIZE_MAX/sizeof(double);
}

static int valid_vector(const struct r64 *array, int64_t length)
{
    return array != NULL && array->data != NULL && array->is_allocated &&
           array->offset >= 0 && length > 0 && length <= INT32_MAX &&
           array->n_dims == 1 && array->dims[0].length == (int32_t)length &&
           array->dims[0].stride == 1;
}

static int valid_matrix3(const struct r64 *array, int64_t columns)
{
    return array != NULL && array->data != NULL && array->is_allocated &&
           array->offset >= 0 && columns > 0 && columns <= INT32_MAX &&
           array->n_dims == 2 && array->dims[0].length == 3 &&
           array->dims[1].length == (int32_t)columns &&
           array->dims[0].stride == (int32_t)columns &&
           array->dims[1].stride == 1;
}

static double *vector_data(struct r64 *array)
{
    return array->data + (size_t)array->offset;
}

static void zero_valid_output(struct r64 *output, int64_t length)
{
    int64_t i;
    if (!valid_vector(output, length)) return;
    for (i = 0; i < length; ++i) vector_data(output)[i] = 0.0;
}

static void interleave_matrix3(const struct r64 *array, int64_t columns,
                               double *interleaved)
{
    const double *data = array->data + (size_t)array->offset;
    int64_t i, j;
    for (j = 0; j < columns; ++j)
        for (i = 0; i < 3; ++i)
            interleaved[3*j+i] =
                data[i*(int64_t)array->dims[0].stride +
                     j*(int64_t)array->dims[1].stride];
}

static int same_matrix3_storage(const struct r64 *left,
                                const struct r64 *right)
{
    return left->data == right->data && left->offset == right->offset &&
           left->n_dims == right->n_dims &&
           left->dims[0].length == right->dims[0].length &&
           left->dims[1].length == right->dims[1].length &&
           left->dims[0].stride == right->dims[0].stride &&
           left->dims[1].stride == right->dims[1].stride;
}

void biesolver_lfmm3d_t_cd_p_rowmajor(
    double eps, int64_t nsource, struct r64 *source,
    struct r64 *charge, struct r64 *dipvec,
    int64_t ntarg, struct r64 *targ, struct r64 *pottarg,
    int64_t *ier, int64_t *adapter_status, double *elapsed_seconds)
{
    size_t source_count, target_count;
    double *source_cm = NULL, *dipvec_cm = NULL, *target_cm = NULL;
    int target_aliases_source = 0;
    double started;

    if (ier != NULL) *ier = 0;
    if (adapter_status != NULL) *adapter_status = 1;
    if (elapsed_seconds != NULL) *elapsed_seconds = 0.0;
    zero_valid_output(pottarg, ntarg);
    if (ier == NULL || adapter_status == NULL || elapsed_seconds == NULL ||
        !(eps > 0.0) || !valid_matrix3(source, nsource) ||
        !valid_vector(charge, nsource) || !valid_matrix3(dipvec, nsource) ||
        !valid_matrix3(targ, ntarg) || !valid_vector(pottarg, ntarg) ||
        !checked_count(nsource, 3, &source_count) ||
        !checked_count(ntarg, 3, &target_count)) return;

    target_aliases_source = nsource == ntarg &&
                            same_matrix3_storage(source, targ);
    started = emscripten_get_now();
    source_cm = (double *)BIESOLVER_FMM3D_MALLOC(
        source_count*sizeof(double));
    dipvec_cm = (double *)BIESOLVER_FMM3D_MALLOC(
        source_count*sizeof(double));
    if (!target_aliases_source)
        target_cm = (double *)BIESOLVER_FMM3D_MALLOC(
            target_count*sizeof(double));
    if (source_cm == NULL || dipvec_cm == NULL ||
        (!target_aliases_source && target_cm == NULL)) goto cleanup;

    interleave_matrix3(source, nsource, source_cm);
    interleave_matrix3(dipvec, nsource, dipvec_cm);
    if (target_aliases_source) {
        target_cm = source_cm;
    } else {
        interleave_matrix3(targ, ntarg, target_cm);
    }
    BIESOLVER_FMM3D_CALL(&eps, &nsource, source_cm, vector_data(charge),
                         dipvec_cm, &ntarg, target_cm,
                         vector_data(pottarg), ier);
    *adapter_status = 0;
    *elapsed_seconds = (emscripten_get_now()-started)/1000.0;

cleanup:
    if (!target_aliases_source) free(target_cm);
    free(dipvec_cm);
    free(source_cm);
}
