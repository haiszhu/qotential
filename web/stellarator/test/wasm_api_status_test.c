#include "wasm_memory_limits.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../native/wasm_api_adapter.c"

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
    const int64_t max_nsrc = (int64_t)
        (BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES / (64 * sizeof(double)));
    assert(cfg->order == 6);
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
    (void)case_result;
    (void)t_fmm_out;
    (void)has_t_fmm;
    (void)t_close_out;
    (void)has_t_close;
    (void)timeinfo_out;
    (void)has_timeinfo;
    *status = biesolver_wasm_memory_preflight(max_nsrc + 1) == 0 ? 15 : 0;
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
    (void)case_result;
}

void biesolver_scope_reset(void) {}

int main(void)
{
    double tgl[8] = {0}, wgl[8] = {0}, dgl[64] = {0};
    double w_bclag[8] = {0}, legmat[64] = {0};
    double umatr[441] = {0}, vmatr[441] = {0};
    assert(solver_simplex_precomp(8, 5, 21, tgl, wgl, dgl, w_bclag,
                                  legmat, umatr, vmatr) == 0);
    assert(solver_prepare_geometry() == 0);
    assert(solver_run_direct() == 0);
    assert(solver_run_close(12, 36, 6, 1, 0.1, tgl, wgl, dgl, w_bclag,
                            legmat, umatr, vmatr) == 215);
    assert(solver_last_error() == 215);
    puts("WASM_STATUS_215_OK");
    return 0;
}
