#ifndef BIESOLVER_WASM_MEMORY_LIMITS_H
#define BIESOLVER_WASM_MEMORY_LIMITS_H

#include <stdint.h>

/* The allocator rejects any single request above this size.  Sharing the
   limit lets the solver refuse an oversized source count before allocating,
   so the caller receives a status instead of an abort. */
#define BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES \
    (UINT64_C(512) * 1024 * 1024)

int biesolver_wasm_memory_preflight(int64_t nsrc);

#endif
