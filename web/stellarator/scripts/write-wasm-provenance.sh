#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PUBLIC_DIR=${PUBLIC_DIR:-$WEB_ROOT/public/wasm}
WASM32_BUILD_DIR=${WASM32_BUILD_DIR:-$WEB_ROOT/build-wasm32}
WASM64_BUILD_DIR=${WASM64_BUILD_DIR:-$WEB_ROOT/build-wasm64}
STAMP32="$WASM32_BUILD_DIR/provenance.stamp"
STAMP64="$WASM64_BUILD_DIR/provenance.stamp"
test -s "$STAMP32" || { echo "missing wasm32 provenance stamp: $STAMP32" >&2; exit 1; }
test -s "$STAMP64" || { echo "missing wasm64 provenance stamp: $STAMP64" >&2; exit 1; }

value() {
  local key=$1 file=$2
  sed -n "s/^${key}=//p" "$file"
}
test "$(value architecture "$STAMP32")" = wasm32
test "$(value architecture "$STAMP64")" = wasm64

common_keys=(
  lfortran_commit lfortran_diff_sha256 lfortran_runtime_patch_sha256
  fmm_fmm3d_commit fmm_manifest_sha256
  fmm_fort2c_commit fmm_fort2c_source_sha256
  fmm_fort2c_patch_sha256 fmm_fort2c_real_patch_sha256
)
for key in "${common_keys[@]}"; do
  left=$(value "$key" "$STAMP32")
  right=$(value "$key" "$STAMP64")
  if [[ -z $left || $left != "$right" ]]; then
    echo "provenance mismatch for $key" >&2
    exit 1
  fi
done

for pair in "solver.js:artifact_js_sha256:$STAMP32" \
            "solver.wasm:artifact_wasm_sha256:$STAMP32" \
            "solver64.js:artifact_js_sha256:$STAMP64" \
            "solver64.wasm:artifact_wasm_sha256:$STAMP64"; do
  IFS=: read -r file key stamp <<<"$pair"
  test -s "$PUBLIC_DIR/$file" || { echo "missing artifact: $file" >&2; exit 1; }
  actual=$(shasum -a 256 "$PUBLIC_DIR/$file" | awk '{print $1}')
  test "$actual" = "$(value "$key" "$stamp")" || {
    echo "artifact hash does not match provenance stamp: $file" >&2
    exit 1
  }
done

next="$PUBLIC_DIR/.PROVENANCE.txt.next.$$"
trap 'rm -f "$next"' EXIT
{
  printf '# Stellarator WebAssembly solver provenance\n\n'
  printf '[wasm32]\n'
  sed 's/[[:space:]]*$//' "$STAMP32"
  printf '\n[wasm64]\n'
  sed 's/[[:space:]]*$//' "$STAMP64"
} > "$next"
mv -f "$next" "$PUBLIC_DIR/PROVENANCE.txt"
trap - EXIT
echo WASM_PROVENANCE_UPDATED
