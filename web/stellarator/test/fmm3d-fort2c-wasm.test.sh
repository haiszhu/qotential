#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
EMCC=${EMCC:-emcc}
BUILD_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-fort2c-wasm.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT

"$WEB_ROOT/scripts/build-fmm3d-wasm.sh" >/dev/null
"$EMCC" -O2 -std=c99 -ffp-contract=off \
  "$SCRIPT_DIR/fmm3d_fort2c_native_test.c" \
  "$WEB_ROOT/build-fmm3d/libfmm3d-wasm.a" -lm \
  -sENVIRONMENT=node -sEXIT_RUNTIME=1 -sSTACK_SIZE=67108864 \
  -sINITIAL_MEMORY=134217728 -sALLOW_MEMORY_GROWTH=1 \
  -o "$BUILD_DIR/fmm3d-fort2c-wasm-test.js"
output=$(node "$BUILD_DIR/fmm3d-fort2c-wasm-test.js")
printf '%s\n' "$output"
grep -q '^FMM3D_FORT2C_NATIVE_CLANG ' <<<"$output"
