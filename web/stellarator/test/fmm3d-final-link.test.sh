#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
: "${LF:?set LF to the patched in-tree lfortran}"
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
LLVM_NM=${LLVM_NM:-llvm-nm}

LF="$LF" FMM3D_SRC="$FMM3D_SRC" "$WEB_ROOT/scripts/build-wasm.sh"

solver_object="$WEB_ROOT/build-wasm/stellarator_solver.o"
adapter_object="$WEB_ROOT/build-wasm/fmm3d_layout_adapter.o"
archive="$WEB_ROOT/build-fmm3d/libfmm3d-wasm.a"
test -s "$solver_object" -a -s "$adapter_object" -a -s "$archive"
undefined=$("$LLVM_NM" -u "$solver_object")
adapter_count=$(grep -Ec '(^|[[:space:]])biesolver_lfmm3d_t_cd_p_rowmajor$' \
  <<<"$undefined" || true)
raw_count=$(grep -Ec '(^|[[:space:]])lfmm3d_t_cd_p_$' <<<"$undefined" || true)
test "$adapter_count" -eq 1
test "$raw_count" -eq 0
test -s "$WEB_ROOT/public/wasm/solver.js"
test -s "$WEB_ROOT/public/wasm/solver.wasm"
echo "FMM3D_FINAL_LINK_OK bytes=$(wc -c < "$WEB_ROOT/public/wasm/solver.wasm" | tr -d ' ')"
