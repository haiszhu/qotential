#include "wasm_memory_limits.h"

#include <assert.h>
#include <stdint.h>
#include <stdio.h>

int main(void)
{
    const int64_t max_nsrc = (int64_t)
        (BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES / (64 * sizeof(double)));
    assert(biesolver_wasm_memory_preflight(max_nsrc) == 1);
    assert(biesolver_wasm_memory_preflight(max_nsrc + 1) == 0);
    assert(biesolver_wasm_memory_preflight(0) == 0);
    puts("WASM_MEMORY_PREFLIGHT_OK");
    return 0;
}
