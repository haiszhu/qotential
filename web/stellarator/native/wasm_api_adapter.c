#include "lfortran_array.h"

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

struct stellarator_case_config {
    int64_t mp, np, order, isimd, ichart;
    bool use_fmm, use_w7x;
    double restol;
};

struct stellarator_case_result {
    int64_t ntri, nsrc, nrender;
    double grf_error;
    struct r64 *sx, *snx, *sw, *ub, *ubn, *u;
    struct r64 *render_xyz, *render_log_error;
    struct i64 *render_triangles;
};

void stellarator_run_case(struct stellarator_case_config *cfg,
                          struct stellarator_case_result *result,
                          int64_t *status,
                          double *t_fmm_out, bool has_t_fmm,
                          double *t_close_out, bool has_t_close,
                          struct r64 *timeinfo_out, bool has_timeinfo);
void stellarator_result_clear(struct stellarator_case_result *result);
void biesolver_scope_reset(void);

enum { EMPTY, GEOMETRY, DIRECT, CLOSE, FINAL };
static int solver_state = EMPTY;
static int solver_error = 0;
static struct r64 sx, snx, sw, ub, ubn, u, render_xyz, render_log_error;
static struct i64 render_triangles;
static struct stellarator_case_result result = {
    0, 0, 0, 0.0, &sx, &snx, &sw, &ub, &ubn, &u,
    &render_xyz, &render_log_error, &render_triangles
};

static ptrdiff_t r64_index(const struct r64 *array, int64_t i, int64_t j)
{
    return (ptrdiff_t)array->offset + (ptrdiff_t)i * array->dims[0].stride +
           (ptrdiff_t)j * array->dims[1].stride;
}

static ptrdiff_t i64_index(const struct i64 *array, int64_t i, int64_t j)
{
    return (ptrdiff_t)array->offset + (ptrdiff_t)i * array->dims[0].stride +
           (ptrdiff_t)j * array->dims[1].stride;
}

static int fail(int code)
{
    solver_error = code;
    return code;
}

int solver_prepare_geometry(void)
{
    if (solver_state != EMPTY) return fail(101);
    solver_error = 0;
    solver_state = GEOMETRY;
    return 0;
}

int solver_run_direct(void)
{
    if (solver_state != GEOMETRY) return fail(102);
    solver_state = DIRECT;
    return 0;
}

int solver_run_close(int64_t mp, int64_t np, int64_t order,
                     int64_t surface, double restol)
{
    struct stellarator_case_config cfg = {mp, np, order, 0, 1, false,
                                          surface == 1, restol};
    int64_t status = 0;
    if (solver_state != DIRECT) return fail(103);
    stellarator_run_case(&cfg, &result, &status,
                         NULL, false, NULL, false, NULL, false);
    if (status != 0) return fail((int)(200 + status));
    solver_state = CLOSE;
    return 0;
}

int solver_finalize_result(void)
{
    if (solver_state != CLOSE) return fail(104);
    solver_state = FINAL;
    return 0;
}

void solver_clear(void)
{
    stellarator_result_clear(&result);
    biesolver_scope_reset();
    solver_state = EMPTY;
    solver_error = 0;
}

int solver_run(int64_t mp, int64_t np, int64_t order,
               int64_t surface, double restol)
{
    int status;
    int64_t hdim;
    uint64_t max_single_allocation_bytes_per_cell;
    solver_clear();
    if (mp <= 0 || np <= 0) return fail(105);
    if (order < 4 || order > 16 || order % 2 != 0) return fail(106);
    if (surface != 0 && surface != 1) return fail(107);
    if (surface == 1 && (!isfinite(restol) || restol <= 0.0)) return fail(108);
    if (surface == 0) {
        /* Built-in surface: at most four triangles per mp*np cell.  The
           largest single allocation is one 64-by-nsrc double-precision
           direct-interaction block.  The W7-X panel count is data-dependent
           and bounded by the chart cap inside the core instead. */
        hdim = order * (order + 1) / 2;
        max_single_allocation_bytes_per_cell =
            UINT64_C(4) * (uint64_t)hdim * 64 * sizeof(double);
        if ((uint64_t)np > (uint64_t)INT64_MAX /
                max_single_allocation_bytes_per_cell ||
            (uint64_t)mp > (uint64_t)INT64_MAX /
                (max_single_allocation_bytes_per_cell * (uint64_t)np))
            return fail(105);
        if ((uint64_t)np > SIZE_MAX / max_single_allocation_bytes_per_cell ||
            (uint64_t)mp > SIZE_MAX /
                (max_single_allocation_bytes_per_cell * (uint64_t)np))
            return fail(105);
    }
    status = solver_prepare_geometry();
    if (status == 0) status = solver_run_direct();
    if (status == 0) status = solver_run_close(mp, np, order, surface, restol);
    if (status == 0) status = solver_finalize_result();
    return status;
}

int64_t solver_result_nsrc(void) { return solver_state == FINAL ? result.nsrc : 0; }
int64_t solver_result_nrender(void) { return solver_state == FINAL ? result.nrender : 0; }
int64_t solver_result_ntriangles(void)
{
    return solver_state == FINAL ? 16 * result.ntri : 0;
}
double solver_result_grf_error(void)
{
    return solver_state == FINAL ? result.grf_error : 0.0;
}
int solver_last_error(void) { return solver_error; }

static int copy_vector(const struct r64 *source, int64_t needed,
                       double *out, int64_t capacity, int error_code)
{
    int64_t i;
    if (solver_state != FINAL || source == NULL || !source->is_allocated)
        return fail(error_code);
    if (out == NULL || capacity < needed) return fail(error_code + 20);
    for (i = 0; i < needed; ++i)
        out[i] = source->data[(ptrdiff_t)source->offset +
                              (ptrdiff_t)i * source->dims[0].stride];
    return 0;
}

static int copy_matrix3(const struct r64 *source, int64_t count,
                        double *out, int64_t capacity, int error_code)
{
    int64_t i, j, needed = 3 * count;
    if (solver_state != FINAL || source == NULL || !source->is_allocated)
        return fail(error_code);
    if (out == NULL || capacity < needed) return fail(error_code + 20);
    for (j = 0; j < count; ++j)
        for (i = 0; i < 3; ++i)
            out[3 * j + i] = source->data[r64_index(source, i, j)];
    return 0;
}

int solver_copy_sx(double *out, int64_t capacity)
{ return copy_matrix3(result.sx, result.nsrc, out, capacity, 301); }
int solver_copy_snx(double *out, int64_t capacity)
{ return copy_matrix3(result.snx, result.nsrc, out, capacity, 302); }
int solver_copy_sw(double *out, int64_t capacity)
{ return copy_vector(result.sw, result.nsrc, out, capacity, 303); }
int solver_copy_ub(double *out, int64_t capacity)
{ return copy_vector(result.ub, result.nsrc, out, capacity, 304); }
int solver_copy_ubn(double *out, int64_t capacity)
{ return copy_vector(result.ubn, result.nsrc, out, capacity, 305); }
int solver_copy_u(double *out, int64_t capacity)
{ return copy_vector(result.u, result.nsrc, out, capacity, 306); }
int solver_copy_render_xyz(double *out, int64_t capacity)
{ return copy_matrix3(result.render_xyz, result.nrender, out, capacity, 307); }
int solver_copy_render_log_error(double *out, int64_t capacity)
{ return copy_vector(result.render_log_error, result.nrender, out, capacity, 308); }

int solver_copy_render_triangles(int64_t *out, int64_t capacity)
{
    int64_t i, j, count = 16 * result.ntri, needed = 3 * count;
    if (solver_state != FINAL || !result.render_triangles->is_allocated)
        return fail(309);
    if (out == NULL || capacity < needed) return fail(329);
    for (j = 0; j < count; ++j)
        for (i = 0; i < 3; ++i)
            out[3 * j + i] = result.render_triangles->data[
                i64_index(result.render_triangles, i, j)];
    return 0;
}
