#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/wasm-artifacts.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT
PUBLIC_DIR="$TMP_DIR/public"
BUILD32="$TMP_DIR/build-wasm32"
BUILD64="$TMP_DIR/build-wasm64"
mkdir -p "$PUBLIC_DIR" "$BUILD32" "$BUILD64"
printf 'old-checksum\n' > "$PUBLIC_DIR/SHA256SUMS"
printf 'old-provenance\n' > "$PUBLIC_DIR/PROVENANCE.txt"
cp "$PUBLIC_DIR/SHA256SUMS" "$TMP_DIR/old.sha"
cp "$PUBLIC_DIR/PROVENANCE.txt" "$TMP_DIR/old.provenance"

for file in solver.js solver.wasm solver64.js; do
  printf '%s\n' "$file" > "$PUBLIC_DIR/$file"
done
: > "$PUBLIC_DIR/solver64.wasm"
if PUBLIC_DIR="$PUBLIC_DIR" "$WEB_ROOT/scripts/update-wasm-checksums.sh"; then
  echo 'checksum writer accepted an empty artifact' >&2
  exit 1
fi
cmp "$TMP_DIR/old.sha" "$PUBLIC_DIR/SHA256SUMS"
printf 'solver64.wasm\n' > "$PUBLIC_DIR/solver64.wasm"

common_stamp() {
  cat <<'EOF'
lfortran_commit=lf-commit
lfortran_diff_sha256=lf-diff
lfortran_runtime_patch_sha256=lf-patch
fmm_fmm3d_commit=fmm-commit
fmm_manifest_sha256=manifest
fmm_fort2c_commit=fort2c-commit
fmm_fort2c_source_sha256=fort2c-source
fmm_fort2c_patch_sha256=fort2c-allocate-patch
fmm_fort2c_real_patch_sha256=fort2c-real-patch
EOF
}
{
  echo architecture=wasm32
  common_stamp
  echo "artifact_js_sha256=$(shasum -a 256 "$PUBLIC_DIR/solver.js" | awk '{print $1}')"
  echo "artifact_wasm_sha256=$(shasum -a 256 "$PUBLIC_DIR/solver.wasm" | awk '{print $1}')"
} > "$BUILD32/provenance.stamp"
{
  echo architecture=wasm64
  common_stamp
  echo fmm_fort2c_source_sha256=mismatch
  echo "artifact_js_sha256=$(shasum -a 256 "$PUBLIC_DIR/solver64.js" | awk '{print $1}')"
  echo "artifact_wasm_sha256=$(shasum -a 256 "$PUBLIC_DIR/solver64.wasm" | awk '{print $1}')"
} > "$BUILD64/provenance.stamp"

if PUBLIC_DIR="$PUBLIC_DIR" WASM32_BUILD_DIR="$BUILD32" WASM64_BUILD_DIR="$BUILD64" \
    "$WEB_ROOT/scripts/write-wasm-provenance.sh"; then
  echo 'provenance writer accepted mismatched translator source' >&2
  exit 1
fi
cmp "$TMP_DIR/old.provenance" "$PUBLIC_DIR/PROVENANCE.txt"
sed -i.bak '/^fmm_fort2c_source_sha256=mismatch$/d' "$BUILD64/provenance.stamp"
rm -f "$BUILD64/provenance.stamp.bak"

PUBLIC_DIR="$PUBLIC_DIR" WASM32_BUILD_DIR="$BUILD32" WASM64_BUILD_DIR="$BUILD64" \
  "$WEB_ROOT/scripts/write-wasm-provenance.sh"
PUBLIC_DIR="$PUBLIC_DIR" "$WEB_ROOT/scripts/update-wasm-checksums.sh"
test "$(wc -l < "$PUBLIC_DIR/SHA256SUMS" | tr -d ' ')" = 4
(cd "$PUBLIC_DIR" && shasum -a 256 --check SHA256SUMS >/dev/null)
grep -q '^\[wasm32\]$' "$PUBLIC_DIR/PROVENANCE.txt"
grep -q '^\[wasm64\]$' "$PUBLIC_DIR/PROVENANCE.txt"
echo WASM_ARTIFACT_POLICY_OK
