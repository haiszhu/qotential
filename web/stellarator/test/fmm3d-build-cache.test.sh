#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
FORT2C=${FORT2C:-fort2c}
FORT2C_SRC=${FORT2C_SRC:-$HOME/fort2c}
EMCC=${EMCC:-emcc}
EMAR=${EMAR:-emar}
LLVM_NM=${LLVM_NM:-llvm-nm}
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-build-cache.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT
BUILDER="$WEB_ROOT/scripts/build-fmm3d-wasm.sh"

run_builder() {
  local memory64=$1
  local build_dir=$2
  local fort2c_src=${3:-$FORT2C_SRC}
  FMM3D_SRC="$FMM3D_SRC" FORT2C="$FORT2C" FORT2C_SRC="$fort2c_src" \
    EMCC="$EMCC" EMAR="$EMAR" LLVM_NM="$LLVM_NM" \
    WASM_MEMORY64="$memory64" FMM3D_BUILD_DIR="$build_dir" \
    "$BUILDER" 2>&1
}

BUILD32="$TMP_DIR/build-fmm3d-wasm32"
first=$(run_builder 0 "$BUILD32")
grep -q 'FMM3D_WASM_REBUILT' <<<"$first"

archive="$BUILD32/libfmm3d-wasm.a"
test -s "$archive"
grep -q '^architecture=wasm32$' "$BUILD32/cache.stamp"
grep -Eq '^fort2c_source_sha256=[0-9a-f]{64}$' "$BUILD32/cache.stamp"
sha_first=$(shasum -a 256 "$archive" | awk '{print $1}')

second=$(run_builder 0 "$BUILD32")
grep -q 'FMM3D_WASM_CACHE_HIT' <<<"$second"
sha_second=$(shasum -a 256 "$archive" | awk '{print $1}')
test "$sha_first" = "$sha_second"

if run_builder 1 "$BUILD32" >"$TMP_DIR/mixed.log" 2>&1; then
  echo 'wasm64 unexpectedly reused the wasm32 cache directory' >&2
  exit 1
fi
grep -q 'belongs to wasm32' "$TMP_DIR/mixed.log"

BUILD64="$TMP_DIR/build-fmm3d-wasm64"
memory64=$(run_builder 1 "$BUILD64")
grep -q 'FMM3D_WASM_REBUILT' <<<"$memory64"
grep -q '^architecture=wasm64$' "$BUILD64/cache.stamp"
grep -q -- '-m64' "$BUILD64/cache.stamp"

FORT2C_COPY="$TMP_DIR/fort2c-copy"
mkdir -p "$FORT2C_COPY"
cp -R "$FORT2C_SRC/fort2c" "$FORT2C_COPY/fort2c"
git -C "$FORT2C_COPY" init -q
git -C "$FORT2C_COPY" add fort2c
git -C "$FORT2C_COPY" -c user.name=cache-test \
  -c user.email=cache-test.invalid commit -q -m baseline
BUILD_MUTATION="$TMP_DIR/build-fort2c-mutation"
mutation_first=$(run_builder 0 "$BUILD_MUTATION" "$FORT2C_COPY")
grep -q 'FMM3D_WASM_REBUILT' <<<"$mutation_first"
mutation_hit=$(run_builder 0 "$BUILD_MUTATION" "$FORT2C_COPY")
grep -q 'FMM3D_WASM_CACHE_HIT' <<<"$mutation_hit"
printf '\n# cache-invalidation probe\n' >> "$FORT2C_COPY/fort2c/__init__.py"
mutation_changed=$(run_builder 0 "$BUILD_MUTATION" "$FORT2C_COPY")
grep -q 'FMM3D_WASM_REBUILT' <<<"$mutation_changed"

members=$("$EMAR" t "$archive")
symbols=$("$LLVM_NM" "$archive")
unresolved=$("$LLVM_NM" -u "$archive")
grep -q 'lfmm3dwrap.o' <<<"$members"
grep -q 'lfmm3d_t_cd_p_' <<<"$symbols"
if grep -q 'lfmm3d_t_cd_p_c_' <<<"$symbols"; then
  echo 'unexpected differential-test FMM symbol in drop-in archive' >&2
  exit 1
fi
if grep -q '_c_$' <<<"$unresolved"; then
  echo 'unresolved differential-test symbol in drop-in archive' >&2
  exit 1
fi

echo "FMM3D_BUILD_CACHE_OK wasm32=$sha_first wasm64=$( \
  shasum -a 256 "$BUILD64/libfmm3d-wasm.a" | awk '{print $1}')"
