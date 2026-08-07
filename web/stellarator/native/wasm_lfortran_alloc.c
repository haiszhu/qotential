#include "wasm_memory_limits.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

/* Generated LFortran C only treats this value as an opaque token.  Bypass
   the runtime's function-pointer allocator table because wasm32 indirect
   function table initialization is not reliable in the pinned compiler. */
struct lfortran_allocator_t;

static uint64_t traced_bytes;
static uint64_t traced_calls;

struct scratch_chunk {
    struct scratch_chunk *next;
    size_t capacity;
    size_t used;
    unsigned char data[];
};

struct biesolver_scope_mark { void *chunk; size_t used; };
void biesolver_scope_reset(void);

static struct scratch_chunk *scratch_first;
static struct scratch_chunk *scratch_current;

void *biesolver_scope_alloc(size_t size)
{
    const size_t alignment = 16;
    size_t aligned = (size + alignment - 1) & ~(alignment - 1);
    size_t capacity;
    struct scratch_chunk *chunk;

    if (aligned < size) return NULL;
    if (scratch_current == NULL ||
        aligned > scratch_current->capacity - scratch_current->used) {
        chunk = scratch_current == NULL ? scratch_first : scratch_current->next;
        if (chunk == NULL || chunk->capacity < aligned) {
            capacity = scratch_current == NULL ? (size_t)1024 * 1024
                                               : scratch_current->capacity * 2;
            if (capacity < aligned) capacity = aligned;
            chunk = malloc(sizeof(*chunk) + capacity);
            if (chunk == NULL) return NULL;
            chunk->next = scratch_current == NULL ? scratch_first
                                                  : scratch_current->next;
            chunk->capacity = capacity;
            chunk->used = 0;
            if (scratch_current == NULL) scratch_first = chunk;
            else scratch_current->next = chunk;
        } else {
            chunk->used = 0;
        }
        scratch_current = chunk;
    }
    chunk = scratch_current;
    void *pointer = chunk->data + chunk->used;
    chunk->used += aligned;
    return pointer;
}

struct biesolver_scope_mark biesolver_scope_enter(void)
{
    return (struct biesolver_scope_mark){
        scratch_current,
        scratch_current == NULL ? 0 : scratch_current->used
    };
}

void biesolver_scope_leave(struct biesolver_scope_mark mark)
{
    struct scratch_chunk *chunk;
    if (mark.chunk == NULL) {
        biesolver_scope_reset();
        return;
    }
    scratch_current = mark.chunk;
    scratch_current->used = mark.used;
    for (chunk = scratch_current->next; chunk != NULL; chunk = chunk->next)
        chunk->used = 0;
}

void biesolver_scope_reset(void)
{
    struct scratch_chunk *chunk;
    for (chunk = scratch_first; chunk != NULL; chunk = chunk->next)
        chunk->used = 0;
    scratch_current = NULL;
}

static void reject_suspicious_allocation(size_t size, const char *origin)
{
    traced_bytes += (uint64_t)size;
    traced_calls += 1;
    /* traced_bytes is cumulative traffic, not live heap usage: this shim
       cannot subtract a size in free(pointer).  Reject only an implausible
       single request; ASan and malloc still diagnose real live-memory bugs. */
    if ((uint64_t)size <= BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES) return;
#ifdef __EMSCRIPTEN__
    emscripten_log(EM_LOG_ERROR | EM_LOG_C_STACK,
                   "suspicious %s allocation(%zu), cumulative=%llu calls=%llu",
                   origin, size, (unsigned long long)traced_bytes,
                   (unsigned long long)traced_calls);
#else
    fprintf(stderr,
            "suspicious %s allocation(%zu), cumulative=%llu calls=%llu\n",
            origin, size, (unsigned long long)traced_bytes,
            (unsigned long long)traced_calls);
#endif
    abort();
}

struct lfortran_allocator_t *_lfortran_get_default_allocator(void)
{
    return NULL;
}

void *_lfortran_malloc_alloc(struct lfortran_allocator_t *unused, int64_t size)
{
    (void)unused;
    if (size < 0) return NULL;
    if ((uint64_t)size > SIZE_MAX) return NULL;
    reject_suspicious_allocation((size_t)size, "LFortran");
    return malloc((size_t)size);
}

void *biesolver_debug_lfortran_alloc(struct lfortran_allocator_t *unused,
                                     int64_t size, int source_line)
{
    (void)unused;
    if (size < 0) return NULL;
    if ((uint64_t)size > SIZE_MAX) return NULL;
    traced_bytes += (uint64_t)size;
    traced_calls += 1;
    if ((uint64_t)size > BIESOLVER_MAX_SINGLE_ALLOCATION_BYTES) {
#ifdef __EMSCRIPTEN__
        emscripten_log(EM_LOG_ERROR,
                       "generated C line %d: LFortran allocation(%lld), cumulative=%llu calls=%llu",
                       source_line, (long long)size,
                       (unsigned long long)traced_bytes,
                       (unsigned long long)traced_calls);
#else
        fprintf(stderr,
                "generated C line %d: LFortran allocation(%lld), cumulative=%llu calls=%llu\n",
                source_line, (long long)size,
                (unsigned long long)traced_bytes,
                (unsigned long long)traced_calls);
#endif
        abort();
    }
    return malloc((size_t)size);
}

void _lfortran_free_alloc(struct lfortran_allocator_t *unused, char *pointer)
{
    (void)unused;
    free(pointer);
}

/* Debug-only entry points used when WASM_ALLOC_TRACE=1.  Only the generated
   solver translation unit has malloc/free remapped to these functions, so the
   implementation itself can safely call the real C allocator. */
void *biesolver_debug_malloc(size_t size)
{
    reject_suspicious_allocation(size, "generated-C");
    return malloc(size);
}

void biesolver_debug_free(void *pointer)
{
    free(pointer);
}
