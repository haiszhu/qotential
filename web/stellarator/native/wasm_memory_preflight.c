/* Reject a source count whose 64-by-nsrc direct-interaction block would
   exceed the allocator's single-request limit.  This is deliberately narrow:
   it prevents the one known abort, and is not a total working-set estimate. */
#include "wasm_memory_limits.h"

int biesolver_wasm_memory_preflight(int64_t nsrc)
{
    const uint64_t bytes_per_source = UINT64_C(64) * sizeof(double);
    if (nsrc <= 0) return 0;
    return (uint64_t)nsrc <=
        BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES / bytes_per_source;
}
