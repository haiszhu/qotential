#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
QOTENTIAL_DIR=$(cd "$WEB_ROOT/../.." && pwd)
LF=${LF:-lfortran}
EMCC=${EMCC:-emcc}
LLVM_NM=${LLVM_NM:-llvm-nm}
WASM_OPT=${WASM_OPT:--O2}
WASM_DEBUG=${WASM_DEBUG:-}
WASM_ALLOC_TRACE=${WASM_ALLOC_TRACE:-0}
LF_ARRAY_BOUNDS_CHECK=${LF_ARRAY_BOUNDS_CHECK:-0}
WASM_SANITIZE=${WASM_SANITIZE:-}
WASM_STACK_SIZE=${WASM_STACK_SIZE:-4194304}
WASM_PROGRESS=${WASM_PROGRESS:-1}
if [[ $WASM_PROGRESS != 0 && $WASM_PROGRESS != 1 ]]; then
  echo "WASM_PROGRESS must be 0 or 1 (got '$WASM_PROGRESS')" >&2
  exit 1
fi
LF_BIN=$(command -v "$LF" 2>/dev/null || printf '%s\n' "$LF")
if [[ ! -x $LF_BIN ]]; then
  echo "LF=$LF is not an executable." >&2
  echo "Build the patched compiler first; see the macOS toolchain section of" >&2
  echo "web/stellarator/README.md.  LF must be <lfortran-source>/build/src/bin/lfortran" >&2
  echo "-- the in-tree binary, because the source tree is recovered from its path." >&2
  exit 1
fi
if [[ -z ${LF_SRC:-} ]]; then
  LF_SRC=$(cd "$(dirname "$LF_BIN")/../../.." && pwd)
fi
RUNTIME_DIR="$LF_SRC/src/libasr/runtime"
if [[ ! -f $RUNTIME_DIR/lfortran_intrinsics.h || ! -f $RUNTIME_DIR/lfortran_intrinsics.c ]]; then
  echo "no LFortran runtime headers under $RUNTIME_DIR" >&2
  echo "LF_SRC was derived as $LF_SRC from LF=$LF_BIN; it must be the LFortran" >&2
  echo "source tree.  Point LF at <source>/build/src/bin/lfortran, or set LF_SRC." >&2
  exit 1
fi

BUILD_DIR="$WEB_ROOT/build-wasm"
MOD_DIR="$BUILD_DIR/modules"
PUBLIC_DIR="$WEB_ROOT/public/wasm"
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR" "$MOD_DIR" "$PUBLIC_DIR"

CPP_FLAGS=(
  --cpp --no-style-suggestions --legacy-array-sections
  -D BIESOLVER_R64_ONLY
  -D BIESOLVER_STELLARATOR_BUILD
  -D BIESOLVER_WASM_SCALAR_ONLY
  -D BIESOLVER_C_BACKEND_ROW_MAJOR
)
if [[ $LF_ARRAY_BOUNDS_CHECK == 0 ]]; then
  CPP_FLAGS+=(--no-array-bounds-checking)
fi
# Progress instrumentation is opt-in and lives entirely behind this macro.  With
# WASM_PROGRESS=0 neither the macro nor the EM_JS bridge object is present, so
# the generated C and solver.wasm reproduce the pre-progress baseline byte-for-byte.
if [[ $WASM_PROGRESS == 1 ]]; then
  CPP_FLAGS+=(-D BIESOLVER_WASM_PROGRESS)
fi
SANITIZER_FLAGS=(-DBIESOLVER_WASM_BUILD=1)
PATH_MAP="-ffile-prefix-map=$WEB_ROOT=."
if [[ -n $WASM_SANITIZE ]]; then
  SANITIZER_FLAGS+=("-fsanitize=$WASM_SANITIZE")
fi

while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  source="$QOTENTIAL_DIR/$relative"
  name=$(basename "${relative%.f90}")
  "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
    -c "$source" -o "$BUILD_DIR/setup-$name.o"
  test -s "$BUILD_DIR/setup-$name.o"
done < "$WEB_ROOT/fortran/wasm_sources.txt"

SOLVER_C="$BUILD_DIR/stellarator_solver.c"
"$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" --show-c \
  "$WEB_ROOT/fortran/stellarator_grf_core_mod.f90" > "$SOLVER_C"
python3 "$WEB_ROOT/scripts/rewrite-lfortran-c.py" "$SOLVER_C"
lines=$(wc -l < "$SOLVER_C" | tr -d ' ')
if (( lines < 1000 )); then
  echo "generated solver C is only $lines lines" >&2
  exit 1
fi

# Not under /tmp: on macOS that is a symlink to /private/tmp, and Emscripten
# builds its cached system libraries with source paths computed relative to the
# cache directory.  The symlink throws that arithmetic off by a level and the
# build dies with "clang: error: no such file or directory:
# '../../../../opt/homebrew/Cellar/emscripten/.../gl.c'".  Any real path works;
# `em-config CACHE` reports Emscripten's own default if you prefer that.
export EM_CACHE=${EM_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/biesolver-emscripten}
SOLVER_CC_FLAGS=(
  -DBIESOLVER_GENERATED_SOLVER=1
  -include "$WEB_ROOT/native/wasm_lfortran_alloc.h"
)
if [[ $WASM_ALLOC_TRACE == 1 ]]; then
  SOLVER_CC_FLAGS+=(
    -Dmalloc=biesolver_debug_malloc
    -Dfree=biesolver_debug_free
  )
fi
"$EMCC" "$WASM_OPT" "$WASM_DEBUG" "$PATH_MAP" -std=c11 -ffp-contract=off -I"$RUNTIME_DIR" \
  "${SANITIZER_FLAGS[@]}" \
  "${SOLVER_CC_FLAGS[@]}" \
  -c "$SOLVER_C" -o "$BUILD_DIR/stellarator_solver.o"
test -s "$BUILD_DIR/stellarator_solver.o"

undef=$($LLVM_NM -u "$BUILD_DIR/stellarator_solver.o")
if grep -Ei 'lq_csimd|csimd|lfmm3d|omp_get_|__kmpc|GOMP_|cpu_time|system_clock' \
     <<<"$undef"; then
  echo "forbidden scalar-build dependency remains" >&2
  exit 1
fi

"$EMCC" "$WASM_OPT" "$WASM_DEBUG" "$PATH_MAP" -std=c11 -D_GNU_SOURCE -I"$LF_SRC/src" -I"$RUNTIME_DIR" \
  "${SANITIZER_FLAGS[@]}" \
  -D_lfortran_get_default_allocator=_lfortran_runtime_get_default_allocator \
  -D_lfortran_malloc_alloc=_lfortran_runtime_malloc_alloc \
  -D_lfortran_free_alloc=_lfortran_runtime_free_alloc \
  -c "$RUNTIME_DIR/lfortran_intrinsics.c" -o "$BUILD_DIR/lfortran_intrinsics.o"
"$EMCC" "$WASM_OPT" "$WASM_DEBUG" "$PATH_MAP" -std=c11 "${SANITIZER_FLAGS[@]}" \
  -c "$WEB_ROOT/native/wasm_lfortran_alloc.c" \
  -o "$BUILD_DIR/wasm_lfortran_alloc.o"

# Under WASM_PROGRESS=1 the EM_JS bridge resolves the Fortran core's
# `biesolver_progress` import at the final link.  The callback stays an internal
# import; it is never added to EXPORTED_FUNCTIONS.
PROGRESS_LINK_OBJS=()
if [[ $WASM_PROGRESS == 1 ]]; then
  "$EMCC" "$WASM_OPT" "$WASM_DEBUG" "$PATH_MAP" -std=c11 "${SANITIZER_FLAGS[@]}" \
    -c "$WEB_ROOT/native/wasm_progress_bridge.c" \
    -o "$BUILD_DIR/wasm_progress_bridge.o"
  test -s "$BUILD_DIR/wasm_progress_bridge.o"
  PROGRESS_LINK_OBJS=("$BUILD_DIR/wasm_progress_bridge.o")
fi

EXPORTED="['_solver_simplex_precomp','_solver_run','_solver_result_nsrc','_solver_result_nrender','_solver_result_ntriangles','_solver_result_grf_error','_solver_copy_sx','_solver_copy_snx','_solver_copy_sw','_solver_copy_ub','_solver_copy_ubn','_solver_copy_u','_solver_copy_render_xyz','_solver_copy_render_log_error','_solver_copy_render_triangles','_solver_last_error','_solver_clear','_malloc','_free']"

"$EMCC" "$WASM_OPT" "$WASM_DEBUG" "$PATH_MAP" -std=c11 -D_GNU_SOURCE -ffp-contract=off \
  "${SANITIZER_FLAGS[@]}" \
  -Wno-unused-variable -I"$LF_SRC/src" -I"$RUNTIME_DIR" -I"$WEB_ROOT/native" \
  "$BUILD_DIR/stellarator_solver.o" \
  "$BUILD_DIR/lfortran_intrinsics.o" \
  "$BUILD_DIR/wasm_lfortran_alloc.o" \
  ${PROGRESS_LINK_OBJS[@]+"${PROGRESS_LINK_OBJS[@]}"} \
  "$QOTENTIAL_DIR/utils/stellarator_geo_mex.c" \
  "$WEB_ROOT/native/geometry_layout_adapter.c" \
  "$WEB_ROOT/native/wasm_blas_shim.c" \
  "$WEB_ROOT/native/wasm_memory_preflight.c" \
  "$WEB_ROOT/native/wasm_api_adapter.c" \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createSolver \
  -sWASM_BIGINT=1 -sALLOW_MEMORY_GROWTH=1 -sENVIRONMENT=web,worker,node \
  -sSTACK_SIZE="$WASM_STACK_SIZE" \
  -sINCOMING_MODULE_JS_API=wasmBinary,locateFile \
  -sEXPORTED_RUNTIME_METHODS=HEAPU8 \
  -sEXPORTED_FUNCTIONS="$EXPORTED" \
  -o "$PUBLIC_DIR/solver.js" -lm

test -s "$PUBLIC_DIR/solver.js"
test -s "$PUBLIC_DIR/solver.wasm"
echo "WASM_BUILD_OK lines=$lines bytes=$(wc -c < "$PUBLIC_DIR/solver.wasm" | tr -d ' ')"
