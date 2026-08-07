#include "wasm_memory_limits.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

#include "../native/wasm_api_adapter.c"

void stellarator_run_case(struct stellarator_case_config *cfg,
                          struct stellarator_case_result *case_result,
                          int64_t *status,
                          double *t_fmm_out, bool has_t_fmm,
                          double *t_close_out, bool has_t_close,
                          struct r64 *timeinfo_out, bool has_timeinfo)
{
    const int64_t max_nsrc = (int64_t)
        (BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES / (64 * sizeof(double)));
    (void)cfg;
    (void)case_result;
    (void)t_fmm_out;
    (void)has_t_fmm;
    (void)t_close_out;
    (void)has_t_close;
    (void)timeinfo_out;
    (void)has_timeinfo;
    *status = biesolver_wasm_memory_preflight(max_nsrc + 1) == 0 ? 15 : 0;
}

void stellarator_result_clear(struct stellarator_case_result *case_result)
{
    (void)case_result;
}

void biesolver_scope_reset(void) {}

int main(void)
{
    assert(solver_prepare_geometry() == 0);
    assert(solver_run_direct() == 0);
    assert(solver_run_close(12, 36, 6, 1, 0.1) == 215);
    assert(solver_last_error() == 215);
    puts("WASM_STATUS_215_OK");
    return 0;
}
