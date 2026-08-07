#ifndef BIESOLVER_WASM_LFORTRAN_ALLOC_H
#define BIESOLVER_WASM_LFORTRAN_ALLOC_H

#include <stddef.h>
#include <stdint.h>
#include <lfortran_intrinsics.h>

void *biesolver_debug_malloc(size_t size);
void biesolver_debug_free(void *pointer);
struct biesolver_scope_mark { void *chunk; size_t used; };
void *biesolver_scope_alloc(size_t size);
struct biesolver_scope_mark biesolver_scope_enter(void);
void biesolver_scope_leave(struct biesolver_scope_mark mark);
void biesolver_scope_reset(void);
void *biesolver_debug_lfortran_alloc(lfortran_allocator_t *allocator,
                                     int64_t size, int source_line);

#ifdef BIESOLVER_WASM_ALLOC_TRACE
#define _lfortran_malloc_alloc(allocator, size) \
    biesolver_debug_lfortran_alloc((allocator), (size), __LINE__)
#endif

#endif
