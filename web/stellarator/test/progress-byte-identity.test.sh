#!/usr/bin/env bash
# Prove that the progress instrumentation is inert when WASM_PROGRESS=0: the
# generated C and the linked wasm must be byte-identical to a baseline taken
# before any progress source existed.
#
#   BASELINE_DIR=<dir> LF=<lfortran> bash test/progress-byte-identity.test.sh capture
#   BASELINE_DIR=<dir> LF=<lfortran> bash test/progress-byte-identity.test.sh verify
#
# The manifest pins the toolchain and every build flag, so a compiler or flag
# change fails before the artifact comparison instead of masquerading as
# leaked progress code.
set -euo pipefail

MODE=${1:-}
case "$MODE" in
  capture|verify) ;;
  *) echo "usage: $0 capture|verify" >&2; exit 2 ;;
esac

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${BASELINE_DIR:?set BASELINE_DIR to the baseline directory}"
: "${LF:?set LF to the patched LFortran executable}"
EMCC=${EMCC:-emcc}

# Every value the build depends on, stated rather than inherited.
export WASM_PROGRESS=0
export WASM_OPT=-O2
export WASM_DEBUG=
export WASM_SANITIZE=
export LF_ARRAY_BOUNDS_CHECK=0
export WASM_ALLOC_TRACE=0
export WASM_STACK_SIZE=4194304

GENERATED_C="$WEB_ROOT/build-wasm/stellarator_solver.c"
WASM_BIN="$WEB_ROOT/public/wasm/solver.wasm"

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/stellarator-byte-identity.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

write_manifest() {
  local out=$1
  {
    echo "LF_PATH=$(command -v "$LF")"
    echo "EMCC_PATH=$(command -v "$EMCC")"
    echo "WASM_PROGRESS=$WASM_PROGRESS"
    echo "WASM_OPT=$WASM_OPT"
    echo "WASM_DEBUG=$WASM_DEBUG"
    echo "WASM_SANITIZE=$WASM_SANITIZE"
    echo "LF_ARRAY_BOUNDS_CHECK=$LF_ARRAY_BOUNDS_CHECK"
    echo "WASM_ALLOC_TRACE=$WASM_ALLOC_TRACE"
    echo "WASM_STACK_SIZE=$WASM_STACK_SIZE"
    echo "--- LF --version ---"
    "$LF" --version
    echo "--- EMCC --version ---"
    "$EMCC" --version
  } > "$out"
}

build() {
  "$WEB_ROOT/scripts/build-wasm.sh" > "$TMP_DIR/build.log" 2>&1 || {
    echo "build failed; see below" >&2
    tail -20 "$TMP_DIR/build.log" >&2
    exit 1
  }
  test -s "$GENERATED_C"
  test -s "$WASM_BIN"
}

if [[ $MODE == capture ]]; then
  if [[ -e $BASELINE_DIR ]]; then
    echo "refusing to overwrite existing baseline: $BASELINE_DIR" >&2
    exit 1
  fi
  build
  mkdir -p "$BASELINE_DIR"
  write_manifest "$BASELINE_DIR/manifest.txt"
  cp "$GENERATED_C" "$BASELINE_DIR/stellarator_solver.c"
  cp "$WASM_BIN" "$BASELINE_DIR/solver.wasm"
  echo "BYTE_BASELINE_CAPTURED $BASELINE_DIR"
  shasum -a 256 "$BASELINE_DIR/stellarator_solver.c" "$BASELINE_DIR/solver.wasm"
  exit 0
fi

test -f "$BASELINE_DIR/manifest.txt" ||
  { echo "no baseline manifest in $BASELINE_DIR" >&2; exit 1; }

# Toolchain first: a flag or compiler mismatch must not look like leaked code.
write_manifest "$TMP_DIR/manifest.txt"
if ! cmp -s "$TMP_DIR/manifest.txt" "$BASELINE_DIR/manifest.txt"; then
  echo "toolchain/flag manifest differs from baseline:" >&2
  diff "$BASELINE_DIR/manifest.txt" "$TMP_DIR/manifest.txt" >&2 || true
  exit 1
fi

build
cmp "$BASELINE_DIR/stellarator_solver.c" "$GENERATED_C" ||
  { echo "generated C differs with WASM_PROGRESS=0" >&2; exit 1; }
cmp "$BASELINE_DIR/solver.wasm" "$WASM_BIN" ||
  { echo "solver.wasm differs with WASM_PROGRESS=0" >&2; exit 1; }
echo "BYTE_IDENTITY_OK"
