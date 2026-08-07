#!/usr/bin/env bash
# Probe the real stellarator dependency closure with LFortran's C backend.
#
# A fresh compiler-only build supplies dependency .mod files without routing
# through host LLVM IR. Success requires exit status zero, non-empty generated
# C, and non-empty host and WebAssembly objects for the final API translation.

set -uo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
QOTENTIAL_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
QA_DIR="$QOTENTIAL_DIR/external/QuatApproximation"
LQ_DIR="$QOTENTIAL_DIR/external/LineQuaaadrature"
LF=${LF:-lfortran}

LF_BIN=$(command -v "$LF" 2>/dev/null || printf '%s\n' "$LF")
if [[ -z ${LF_SRC:-} ]]; then
  LF_SRC=$(cd "$(dirname "$LF_BIN")/../../.." 2>/dev/null && pwd)
fi
LF_RUNTIME_DIR="$LF_SRC/src/libasr/runtime"
if [[ ! -f "$LF_RUNTIME_DIR/lfortran_intrinsics.h" ]]; then
  echo "SETUP FAIL: set LF_SRC to the LFortran source tree containing src/libasr/runtime"
  exit 2
fi

PROBE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/qotential-cprobe.XXXXXX")
trap 'rm -rf "$PROBE_DIR"' EXIT

MOD_DIR="$PROBE_DIR/modules"
mkdir -p "$MOD_DIR"
CPP_FLAGS=(
  --cpp
  --no-style-suggestions
  --legacy-array-sections
  --no-array-bounds-checking
  -D BIESOLVER_R64_ONLY
  -D BIESOLVER_STELLARATOR_BUILD
)

NAMES=(
  quatapproximation harmonic koorn_geom qkernel omega tensor_geom
  linequaaadrature adaptive kernel solidangle lap3d lap3d_simd
  patch_refine lap3d_close stellarator_wasm_api
)
SOURCES=(
  "$QA_DIR/src/quatapproximation_mod.f90"
  "$QA_DIR/src/harmonic_mod.f90"
  "$QA_DIR/src/koorn_geom_mod.f90"
  "$QA_DIR/src/qkernel_mod.f90"
  "$QA_DIR/src/omega_mod.f90"
  "$QA_DIR/src/tensor_geom_mod.f90"
  "$LQ_DIR/src/linequaaadrature_mod.f90"
  "$LQ_DIR/src/adaptive_mod.f90"
  "$LQ_DIR/src/kernel_mod.f90"
  "$LQ_DIR/src/solidangle_mod.f90"
  "$LQ_DIR/src/lap3d_mod.f90"
  "$LQ_DIR/src/lap3d_simd_mod.f90"
  "$QOTENTIAL_DIR/src/patch_refine_mod.f90"
  "$QOTENTIAL_DIR/src/lap3d_close_mod.f90"
  "$SCRIPT_DIR/../fortran/stellarator_wasm_api.f90"
)

for i in "${!NAMES[@]}"; do
  name=${NAMES[$i]}
  source=${SOURCES[$i]}
  obj="$PROBE_DIR/mod-$name.o"
  err="$PROBE_DIR/mod-$name.err"
  if ! "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
      -c "$source" -o "$obj" >"$PROBE_DIR/mod-$name.out" 2>"$err"; then
    echo "SETUP FAIL: fresh module build failed for $name"
    tail -20 "$err"
    exit 2
  fi
  if [[ ! -s $obj ]]; then
    echo "SETUP FAIL: module build produced an empty object for $name"
    exit 2
  fi
done

fail=0
printf '%-20s %5s %7s  %s\n' module rc c_lines diagnostic

for i in "${!NAMES[@]}"; do
  name=${NAMES[$i]}
  source=${SOURCES[$i]}
  c_file="$PROBE_DIR/$name.c"
  err_file="$PROBE_DIR/$name.err"

  "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
    --show-c "$source" >"$c_file" 2>"$err_file"
  rc=$?
  lines=$(wc -l <"$c_file" | tr -d ' ')
  diagnostic=$(sed 's/\x1b\[[0-9;]*m//g' "$err_file" \
    | grep -E 'LCompilersException:|code generation error|not supported' \
    | tail -1)

  if [[ $rc -ne 0 || $lines -lt 20 ]]; then
    fail=1
    if [[ -z $diagnostic ]]; then
      diagnostic='<no diagnostic>'
    fi
  else
    diagnostic='OK'
  fi

  printf '%-20s %5d %7d  %s\n' "$name" "$rc" "$lines" "$diagnostic"
done

if [[ $fail -eq 0 ]]; then
  api_c="$PROBE_DIR/stellarator_wasm_api.c"
  if ! clang -std=c11 -I"$LF_RUNTIME_DIR" -c "$api_c" \
      -o "$PROBE_DIR/stellarator-api-host.o" \
      >"$PROBE_DIR/clang.out" 2>"$PROBE_DIR/clang.err"; then
    echo "HOST OBJECT FAIL"
    tail -20 "$PROBE_DIR/clang.err"
    fail=1
  elif [[ ! -s "$PROBE_DIR/stellarator-api-host.o" ]]; then
    echo "HOST OBJECT FAIL: empty object"
    fail=1
  fi

  if ! emcc -std=c11 -I"$LF_RUNTIME_DIR" -c "$api_c" \
      -o "$PROBE_DIR/stellarator-api-wasm.o" \
      >"$PROBE_DIR/emcc.out" 2>"$PROBE_DIR/emcc.err"; then
    echo "WASM OBJECT FAIL"
    tail -20 "$PROBE_DIR/emcc.err"
    fail=1
  elif [[ ! -s "$PROBE_DIR/stellarator-api-wasm.o" ]]; then
    echo "WASM OBJECT FAIL: empty object"
    fail=1
  fi

  if [[ $fail -eq 0 ]]; then
    echo "C_AND_WASM_OBJECT_GATE_OK"
  fi
fi

exit "$fail"
