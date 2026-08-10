#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
EMCC=${EMCC:-emcc}
WASM_MEMORY64=${WASM_MEMORY64:-0}
if [[ $WASM_MEMORY64 == 1 ]]; then
  COMPILE_FLAGS=(-O2 -m64)
  MEMORY_FLAGS=(-sINITIAL_MEMORY=134217728 -sMAXIMUM_MEMORY=17179869184)
  FMM3D_BUILD_DIR=${FMM3D_BUILD_DIR:-$WEB_ROOT/build-fmm3d-wasm64}
else
  COMPILE_FLAGS=(-O2)
  MEMORY_FLAGS=(-sINITIAL_MEMORY=134217728)
  FMM3D_BUILD_DIR=${FMM3D_BUILD_DIR:-$WEB_ROOT/build-fmm3d-wasm32}
fi
BUILD_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-fort2c-wasm.XXXXXX")
TEST_EM_CACHE="$WEB_ROOT/build-emscripten-cache-test.$$"
trap 'rm -rf "$BUILD_DIR" "$TEST_EM_CACHE"' EXIT
export EM_CACHE=${EM_CACHE:-$TEST_EM_CACHE}
mkdir -p "$EM_CACHE"

WASM_MEMORY64="$WASM_MEMORY64" FMM3D_BUILD_DIR="$FMM3D_BUILD_DIR" \
  "$WEB_ROOT/scripts/build-fmm3d-wasm.sh" >/dev/null
"$EMCC" "${COMPILE_FLAGS[@]}" -std=c99 -ffp-contract=off \
  "$SCRIPT_DIR/fmm3d_fort2c_native_test.c" \
  "$FMM3D_BUILD_DIR/libfmm3d-wasm.a" -lm \
  -sENVIRONMENT=node -sEXIT_RUNTIME=1 -sSTACK_SIZE=67108864 \
  "${MEMORY_FLAGS[@]}" -sALLOW_MEMORY_GROWTH=1 \
  -o "$BUILD_DIR/fmm3d-fort2c-wasm-test.js"
output=$(node "$BUILD_DIR/fmm3d-fort2c-wasm-test.js")
printf '%s\n' "$output"
grep -q '^FMM3D_FORT2C_NATIVE_CLANG ' <<<"$output"
