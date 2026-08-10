#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${LF:?set LF to the patched in-tree lfortran}"
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
LLVM_NM=${LLVM_NM:-llvm-nm}
WASM_MEMORY64=${WASM_MEMORY64:-0}
if [[ $WASM_MEMORY64 == 1 ]]; then
  build_dir="$WEB_ROOT/build-wasm64"
  fmm_dir="$WEB_ROOT/build-fmm3d-wasm64"
  basename=solver64
else
  build_dir="$WEB_ROOT/build-wasm32"
  fmm_dir="$WEB_ROOT/build-fmm3d-wasm32"
  basename=solver
fi

LF="$LF" FMM3D_SRC="$FMM3D_SRC" WASM_MEMORY64="$WASM_MEMORY64" \
  "$WEB_ROOT/scripts/build-wasm.sh"

solver_object="$build_dir/stellarator_solver.o"
adapter_object="$build_dir/fmm3d_layout_adapter.o"
archive="$fmm_dir/libfmm3d-wasm.a"
test -s "$solver_object" -a -s "$adapter_object" -a -s "$archive"
undefined=$("$LLVM_NM" -u "$solver_object")
adapter_count=$(grep -Ec '(^|[[:space:]])biesolver_lfmm3d_t_cd_p_rowmajor$' \
  <<<"$undefined" || true)
raw_count=$(grep -Ec '(^|[[:space:]])lfmm3d_t_cd_p_$' <<<"$undefined" || true)
test "$adapter_count" -eq 1
test "$raw_count" -eq 0
test -s "$WEB_ROOT/public/wasm/$basename.js"
test -s "$WEB_ROOT/public/wasm/$basename.wasm"
echo "FMM3D_FINAL_LINK_OK architecture=$WASM_MEMORY64 bytes=$(wc -c < "$WEB_ROOT/public/wasm/$basename.wasm" | tr -d ' ')"
