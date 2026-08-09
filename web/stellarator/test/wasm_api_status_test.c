#include "wasm_memory_limits.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../native/wasm_api_adapter.c"

static int64_t next_core_status;
static int core_call_count;
static int clear_call_count;
static bool expected_use_fmm;
static double expected_fmm_eps;

void stellarator_run_case(struct stellarator_case_config *cfg,
                          struct r64 *tgl, struct r64 *wgl,
                          struct r64 *dgl, struct r64 *w_bclag,
                          struct r64 *legmat, struct r64 *umatr,
                          struct r64 *vmatr,
                          struct stellarator_case_result *case_result,
                          int64_t *status,
                          double *t_fmm_out, bool has_t_fmm,
                          double *t_close_out, bool has_t_close,
                          struct r64 *timeinfo_out, bool has_timeinfo)
{
    ++core_call_count;
    assert(cfg->order == 6);
    assert(cfg->use_fmm == expected_use_fmm);
    assert(cfg->fmm_eps == expected_fmm_eps);
    assert(tgl->n_dims == 1 && tgl->dims[0].length == 8 &&
           tgl->dims[0].stride == 1 && tgl->data[0] == 11.0);
    assert(wgl->n_dims == 1 && wgl->data[0] == 12.0);
    assert(dgl->n_dims == 2 && dgl->dims[0].length == 8 &&
           dgl->dims[1].length == 8 && dgl->dims[0].stride == 8 &&
           dgl->dims[1].stride == 1 && dgl->data[0] == 13.0);
    assert(w_bclag->n_dims == 1 && w_bclag->data[0] == 14.0);
    assert(legmat->n_dims == 2 && legmat->data[0] == 15.0);
    assert(umatr->n_dims == 2 && umatr->dims[0].length == 21 &&
           umatr->dims[1].length == 21 && umatr->data[0] == 16.0);
    assert(vmatr->n_dims == 2 && vmatr->data[0] == 17.0);
    case_result->sx->is_allocated = true;
    (void)t_fmm_out;
    (void)has_t_fmm;
    (void)t_close_out;
    (void)has_t_close;
    (void)timeinfo_out;
    (void)has_timeinfo;
    *status = next_core_status;
}

void simplex_precomp_r64(int64_t nquad, int64_t korder, int64_t kpols,
                         struct r64 *tgl, struct r64 *wgl, struct r64 *dgl,
                         struct r64 *w_bclag, struct r64 *legmat,
                         struct r64 *umatr, struct r64 *vmatr)
{
    assert(nquad == 8 && korder == 5 && kpols == 21);
    assert(tgl->dims[0].length == 8 && tgl->dims[0].stride == 1);
    assert(dgl->dims[0].stride == 8 && dgl->dims[1].stride == 1);
    assert(umatr->dims[0].length == 21 && umatr->dims[1].length == 21);
    tgl->data[0] = 11.0;
    wgl->data[0] = 12.0;
    dgl->data[0] = 13.0;
    w_bclag->data[0] = 14.0;
    legmat->data[0] = 15.0;
    umatr->data[0] = 16.0;
    vmatr->data[0] = 17.0;
}

void stellarator_result_clear(struct stellarator_case_result *case_result)
{
    ++clear_call_count;
    case_result->sx->is_allocated = false;
}

void biesolver_scope_reset(void) {}

int main(void)
{
    static const double tolerances[] = {1e-3, 1e-6, 1e-9, 1e-12, 1e-15};
    double tgl[8] = {0}, wgl[8] = {0}, dgl[64] = {0};
    double w_bclag[8] = {0}, legmat[64] = {0};
    double umatr[441] = {0}, vmatr[441] = {0};
    assert(solver_simplex_precomp(8, 5, 21, tgl, wgl, dgl, w_bclag,
                                  legmat, umatr, vmatr) == 0);
    int i, kernel, before;

    before = core_call_count;
    assert(solver_run(12, 36, 6, 0, 0.1, -1, 1e-3,
        tgl, wgl, dgl, w_bclag, legmat, umatr, vmatr) == 110);
    assert(solver_run(12, 36, 6, 0, 0.1, 2, 1e-3,
        tgl, wgl, dgl, w_bclag, legmat, umatr, vmatr) == 110);
    assert(solver_run(12, 36, 6, 0, 0.1, 1, 2e-4,
        tgl, wgl, dgl, w_bclag, legmat, umatr, vmatr) == 111);
    assert(core_call_count == before);

    next_core_status = 0;
    for (kernel = 0; kernel <= 1; ++kernel) {
        expected_use_fmm = kernel == 1;
        for (i = 0; i < 5; ++i) {
            expected_fmm_eps = tolerances[i];
            assert(solver_run(12, 36, 6, 0, 0.1, kernel, tolerances[i],
                tgl, wgl, dgl, w_bclag, legmat, umatr, vmatr) == 0);
        }
    }

    before = core_call_count;
    assert(solver_run(INT64_C(1000000000), INT64_C(1000000000),
        6, 0, 0.1, 0, 1e-3, tgl, wgl, dgl, w_bclag,
        legmat, umatr, vmatr) == 105);
    assert(core_call_count == before);
    expected_use_fmm = true;
    expected_fmm_eps = 1e-3;
    assert(solver_run(INT64_C(1000000000), INT64_C(1000000000),
        6, 0, 0.1, 1, 1e-3, tgl, wgl, dgl, w_bclag,
        legmat, umatr, vmatr) == 0);
    assert(core_call_count == before + 1);

    for (i = 0; i < 4; ++i) {
        static const int64_t core_status[] = {16, 17, 18, 19};
        int clear_before = clear_call_count;
        next_core_status = core_status[i];
        expected_use_fmm = true;
        expected_fmm_eps = 1e-6;
        assert(solver_run(12, 36, 6, 0, 0.1, 1, 1e-6,
            tgl, wgl, dgl, w_bclag, legmat, umatr, vmatr) == 216 + i);
        assert(solver_last_error() == 216 + i);
        assert(clear_call_count == clear_before + 2);
        assert(!result.sx->is_allocated);
    }
    puts("WASM_FMM_API_STATUS_OK");
    return 0;
}
