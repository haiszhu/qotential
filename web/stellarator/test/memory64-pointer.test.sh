#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
EMCC=${EMCC:-emcc}
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/memory64-pointer.XXXXXX")
TEST_EM_CACHE="$WEB_ROOT/build-emscripten-cache-test.$$"
trap 'rm -rf "$TMP_DIR" "$TEST_EM_CACHE"' EXIT
export EM_CACHE=${EM_CACHE:-$TEST_EM_CACHE}
mkdir -p "$EM_CACHE"

"$EMCC" -O2 -m64 "$WEB_ROOT/test/memory64_pointer_probe.c" \
  -sWASM_BIGINT=1 -sALLOW_MEMORY_GROWTH=1 \
  -sINITIAL_MEMORY=134217728 -sMAXIMUM_MEMORY=17179869184 \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createMemory64Probe \
  -sENVIRONMENT=node -sINCOMING_MODULE_JS_API=wasmBinary,locateFile \
  -sEXPORTED_RUNTIME_METHODS=HEAPU8 \
  -sEXPORTED_FUNCTIONS="['_pointer_identity','_pointer_store','_malloc','_free']" \
  -o "$TMP_DIR/pointer-probe.mjs"

test -s "$TMP_DIR/pointer-probe.mjs"
test -s "$TMP_DIR/pointer-probe.wasm"
node "$WEB_ROOT/test/memory64-pointer.runner.mjs" "$TMP_DIR/pointer-probe.mjs"
