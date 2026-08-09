#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
FORT2C=${FORT2C:-fort2c}
EMCC=${EMCC:-emcc}
EMAR=${EMAR:-emar}
LLVM_NM=${LLVM_NM:-llvm-nm}
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-build-cache.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT
BUILD_DIR="$TMP_DIR/build-fmm3d"
BUILDER="$WEB_ROOT/scripts/build-fmm3d-wasm.sh"

first=$(FMM3D_SRC="$FMM3D_SRC" FORT2C="$FORT2C" EMCC="$EMCC" \
  EMAR="$EMAR" FMM3D_BUILD_DIR="$BUILD_DIR" "$BUILDER" 2>&1)
grep -q 'FMM3D_WASM_REBUILT' <<<"$first"

archive="$BUILD_DIR/libfmm3d-wasm.a"
test -s "$archive"
sha_first=$(shasum -a 256 "$archive" | awk '{print $1}')

second=$(FMM3D_SRC="$FMM3D_SRC" FORT2C="$FORT2C" EMCC="$EMCC" \
  EMAR="$EMAR" FMM3D_BUILD_DIR="$BUILD_DIR" "$BUILDER" 2>&1)
grep -q 'FMM3D_WASM_CACHE_HIT' <<<"$second"
sha_second=$(shasum -a 256 "$archive" | awk '{print $1}')
test "$sha_first" = "$sha_second"

forced=$(FMM3D_SRC="$FMM3D_SRC" FORT2C="$FORT2C" EMCC="$EMCC" \
  EMAR="$EMAR" FMM3D_BUILD_DIR="$BUILD_DIR" FMM3D_REBUILD=1 \
  "$BUILDER" 2>&1)
grep -q 'FMM3D_WASM_REBUILT' <<<"$forced"

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

echo "FMM3D_BUILD_CACHE_OK sha256=$sha_first"
