#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
EMCC=${EMCC:-emcc}
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
WASM_MEMORY64=${WASM_MEMORY64:-0}
if [[ $WASM_MEMORY64 == 1 ]]; then
  COMPILE_FLAGS=(-O2 -m64)
  ALLOC_FLAGS=(-O0 -m64)
  MEMORY_FLAGS=(-sALLOW_MEMORY_GROWTH=1 -sINITIAL_MEMORY=134217728
    -sMAXIMUM_MEMORY=17179869184)
  BUILD_DIR=${FMM3D_BUILD_DIR:-$WEB_ROOT/build-fmm3d-wasm64}
else
  COMPILE_FLAGS=(-O2)
  ALLOC_FLAGS=(-O0)
  MEMORY_FLAGS=(-sALLOW_MEMORY_GROWTH=1)
  BUILD_DIR=${FMM3D_BUILD_DIR:-$WEB_ROOT/build-fmm3d-wasm32}
fi
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-layout-adapter.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT
export EM_CACHE=${EM_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/biesolver-emscripten}

FMM3D_SRC="$FMM3D_SRC" FMM3D_BUILD_DIR="$BUILD_DIR" \
  WASM_MEMORY64="$WASM_MEMORY64" \
  "$WEB_ROOT/scripts/build-fmm3d-wasm.sh" >/dev/null

"$EMCC" "${COMPILE_FLAGS[@]}" -std=c11 -ffp-contract=off \
  -I"$WEB_ROOT/native" \
  "$WEB_ROOT/test/fmm3d_layout_adapter_test.c" \
  "$WEB_ROOT/native/fmm3d_layout_adapter.c" \
  "$BUILD_DIR/libfmm3d-wasm.a" \
  -sSTACK_SIZE=67108864 "${MEMORY_FLAGS[@]}" -sEXIT_RUNTIME=1 \
  -o "$TMP_DIR/layout-parity.js" -lm
node "$TMP_DIR/layout-parity.js"

"$EMCC" "${ALLOC_FLAGS[@]}" -std=c11 -ffp-contract=off \
  -DTEST_ALLOCATION_FAILURE_ONLY \
  -DBIESOLVER_FMM3D_MALLOC=biesolver_test_malloc \
  -DBIESOLVER_FMM3D_CALL=biesolver_test_fmm_call \
  -I"$WEB_ROOT/native" \
  "$WEB_ROOT/test/fmm3d_layout_adapter_test.c" \
  "$WEB_ROOT/native/fmm3d_layout_adapter.c" \
  "${MEMORY_FLAGS[@]}" \
  -sEXIT_RUNTIME=1 -o "$TMP_DIR/layout-allocation.js" -lm
node "$TMP_DIR/layout-allocation.js"
