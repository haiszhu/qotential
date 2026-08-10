#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)

FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
FORT2C=${FORT2C:-fort2c}
FORT2C_SRC=${FORT2C_SRC:-$HOME/fort2c}
EMCC=${EMCC:-emcc}
EMAR=${EMAR:-emar}
LLVM_NM=${LLVM_NM:-llvm-nm}
WASM_MEMORY64=${WASM_MEMORY64:-0}
if [[ $WASM_MEMORY64 != 0 && $WASM_MEMORY64 != 1 ]]; then
  echo "WASM_MEMORY64 must be 0 or 1 (got '$WASM_MEMORY64')" >&2
  exit 2
fi
if [[ $WASM_MEMORY64 == 1 ]]; then
  architecture=wasm64
  default_build_dir="$WEB_ROOT/build-fmm3d-wasm64"
else
  architecture=wasm32
  default_build_dir="$WEB_ROOT/build-fmm3d-wasm32"
fi
FMM3D_BUILD_DIR=${FMM3D_BUILD_DIR:-$default_build_dir}
FMM3D_REBUILD=${FMM3D_REBUILD:-0}
FMM3D_CFLAGS_EXTRA=${FMM3D_CFLAGS_EXTRA:-}

MANIFEST="$WEB_ROOT/fortran/fmm3d_wasm_sources.txt"
RUNTIME_HEADER="$WEB_ROOT/native/fmm3d_c.h"
FORT2C_PATCH="$WEB_ROOT/patches/fort2c-allocate-stat.patch"
FORT2C_REAL_PATCH="$WEB_ROOT/patches/fort2c-real-complex-kind.patch"
FORT2C_RUNNER="$WEB_ROOT/scripts/run-fort2c-source.sh"
SOURCE_DIGEST="$WEB_ROOT/scripts/source-tree-digest.py"
ARCHIVE="$FMM3D_BUILD_DIR/libfmm3d-wasm.a"
STAMP="$FMM3D_BUILD_DIR/cache.stamp"

sha256_file() {
  shasum -a 256 "$1" | awk '{print $1}'
}

require_tool() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "required tool not found: $1" >&2
    exit 1
  fi
}

require_tool "$FORT2C"
require_tool "$EMCC"
require_tool "$EMAR"
require_tool "$LLVM_NM"
require_tool git
require_tool shasum
require_tool python3

test -f "$MANIFEST" || { echo "missing FMM3D manifest: $MANIFEST" >&2; exit 1; }
test -f "$RUNTIME_HEADER" || { echo "missing fort2c runtime header: $RUNTIME_HEADER" >&2; exit 1; }
test -f "$FORT2C_PATCH" || { echo "missing fort2c patch: $FORT2C_PATCH" >&2; exit 1; }
test -f "$FORT2C_REAL_PATCH" || { echo "missing fort2c patch: $FORT2C_REAL_PATCH" >&2; exit 1; }
test -x "$FORT2C_RUNNER" || { echo "missing fort2c source runner: $FORT2C_RUNNER" >&2; exit 1; }
test -x "$SOURCE_DIGEST" || { echo "missing source digest helper: $SOURCE_DIGEST" >&2; exit 1; }
test -d "$FMM3D_SRC/.git" || { echo "FMM3D_SRC is not a git checkout: $FMM3D_SRC" >&2; exit 1; }
test -d "$FORT2C_SRC/.git" || { echo "FORT2C_SRC is not a git checkout: $FORT2C_SRC" >&2; exit 1; }
test -f "$FORT2C_SRC/fort2c/transpiler.py" || {
  echo "missing fort2c package under $FORT2C_SRC" >&2
  exit 1
}

fmm_tag=$(git -C "$FMM3D_SRC" describe --tags --match 'v[0-9]*' --abbrev=0 2>/dev/null || true)
version=${fmm_tag#v}
major=${version%%.*}
minor_tail=${version#*.}
minor=${minor_tail%%.*}
if [[ -z $fmm_tag || ! $major =~ ^[0-9]+$ || ! $minor =~ ^[0-9]+$ ]] || \
   (( major < 2 )); then
  echo "FMM3D >= 2.0 is required; nearest version tag is '${fmm_tag:-none}'" >&2
  exit 1
fi

FMM3D_SRC="$FMM3D_SRC" bash "$WEB_ROOT/test/fmm3d-manifest.test.sh" >/dev/null
FORT2C="$FORT2C_RUNNER" FORT2C_TOOL="$FORT2C" FORT2C_SRC="$FORT2C_SRC" \
  bash "$WEB_ROOT/test/fort2c-allocate-stat.test.sh" >/dev/null
FORT2C="$FORT2C_RUNNER" FORT2C_TOOL="$FORT2C" FORT2C_SRC="$FORT2C_SRC" \
  FORT2C_CHECK_BASELINE=0 \
  bash "$WEB_ROOT/test/fort2c-real-complex-kind.test.sh" >/dev/null

extra_flags=()
if [[ -n $FMM3D_CFLAGS_EXTRA ]]; then
  read -r -a extra_flags <<<"$FMM3D_CFLAGS_EXTRA"
fi
base_flags=(
  -O2 -std=c99 -DFMM3D_DROP_IN -ffp-contract=off
  -fno-strict-aliasing -fwrapv
  -Wno-unused-parameter -Wno-unused-variable -Wno-parentheses
  -Wno-deprecated-non-prototype
)
if [[ $WASM_MEMORY64 == 1 ]]; then
  base_flags+=(-m64)
fi

fmm_commit=$(git -C "$FMM3D_SRC" rev-parse HEAD)
fmm_describe=$(git -C "$FMM3D_SRC" describe --tags --always --dirty)
fort2c_commit=$(git -C "$FORT2C_SRC" rev-parse HEAD)
fort2c_describe=$(git -C "$FORT2C_SRC" describe --tags --always --dirty)
fort2c_path=$(command -v "$FORT2C")
fort2c_version=$(FORT2C_TOOL="$FORT2C" FORT2C_SRC="$FORT2C_SRC" \
  "$FORT2C_RUNNER" --version)
fort2c_source_sha=$(python3 "$SOURCE_DIGEST" "$FORT2C_SRC/fort2c")
emcc_version=$("$EMCC" --version | sed -n '1p')
manifest_sha=$(sha256_file "$MANIFEST")
header_sha=$(sha256_file "$RUNTIME_HEADER")
patch_sha=$(sha256_file "$FORT2C_PATCH")
real_patch_sha=$(sha256_file "$FORT2C_REAL_PATCH")

source_digests=""
while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  source="$FMM3D_SRC/$relative"
  test -f "$source" || { echo "missing FMM3D source: $source" >&2; exit 1; }
  source_digests+="$relative=$(sha256_file "$source")"$'\n'
done < "$MANIFEST"

mkdir -p "$FMM3D_BUILD_DIR"
if [[ -f $STAMP ]]; then
  cached_architecture=$(sed -n 's/^architecture=//p' "$STAMP")
  if [[ -z $cached_architecture ]]; then
    echo "existing FMM3D cache stamp has no architecture: $STAMP" >&2
    exit 1
  fi
  if [[ $cached_architecture != "$architecture" ]]; then
    echo "FMM3D cache directory belongs to $cached_architecture, not $architecture: $FMM3D_BUILD_DIR" >&2
    exit 1
  fi
fi
work_dir=$(mktemp -d "$FMM3D_BUILD_DIR/.build.XXXXXX")
trap 'rm -rf "$work_dir"' EXIT
gen_dir="$work_dir/generated"
obj_dir="$work_dir/objects"
mkdir -p "$gen_dir" "$obj_dir"
cp "$RUNTIME_HEADER" "$gen_dir/fmm3d_c.h"

new_stamp="$work_dir/cache.stamp"
{
  printf 'architecture=%s\n' "$architecture"
  printf 'fmm3d_commit=%s\n' "$fmm_commit"
  printf 'fmm3d_describe=%s\n' "$fmm_describe"
  printf 'fort2c_path=%s\n' "$fort2c_path"
  printf 'fort2c_version=%s\n' "$fort2c_version"
  printf 'fort2c_commit=%s\n' "$fort2c_commit"
  printf 'fort2c_describe=%s\n' "$fort2c_describe"
  printf 'fort2c_source_sha256=%s\n' "$fort2c_source_sha"
  printf 'fort2c_patch_sha256=%s\n' "$patch_sha"
  printf 'fort2c_real_patch_sha256=%s\n' "$real_patch_sha"
  printf 'emcc_version=%s\n' "$emcc_version"
  printf 'manifest_sha256=%s\n' "$manifest_sha"
  printf 'runtime_header_sha256=%s\n' "$header_sha"
  printf 'cflags='
  separator=
  for flag in "${base_flags[@]}"; do
    printf '%s%q' "$separator" "$flag"
    separator=' '
  done
  if [[ -n $FMM3D_CFLAGS_EXTRA ]]; then
    for flag in "${extra_flags[@]}"; do
      printf '%s%q' "$separator" "$flag"
      separator=' '
    done
  fi
  printf '\n'
  printf '%s' "$source_digests"
} > "$new_stamp"

if [[ $FMM3D_REBUILD != 1 && -s $ARCHIVE && -f $STAMP ]] && \
   cmp -s "$new_stamp" "$STAMP"; then
  echo "FMM3D_WASM_CACHE_HIT archive=$ARCHIVE commit=$fmm_commit"
  exit 0
fi

objects=()
while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  source="$FMM3D_SRC/$relative"
  name=$(basename "${relative%.f}")
  generated="$gen_dir/$name.c"
  object="$obj_dir/$name.o"
  FORT2C_TOOL="$FORT2C" FORT2C_SRC="$FORT2C_SRC" \
    "$FORT2C_RUNNER" "$source" --runtime-header fmm3d_c.h \
    --guard-prefix BIESOLVER_FMM3D_ -o "$gen_dir" >/dev/null
  test -s "$generated" || { echo "fort2c produced no C for $relative" >&2; exit 1; }
  "$EMCC" "${base_flags[@]}" \
    ${extra_flags[@]+"${extra_flags[@]}"} -I"$gen_dir" \
    -c "$generated" -o "$object"
  test -s "$object" || { echo "emcc produced no object for $relative" >&2; exit 1; }
  objects+=("$object")
done < "$MANIFEST"

new_archive="$work_dir/libfmm3d-wasm.a"
"$EMAR" rcs "$new_archive" "${objects[@]}"
test -s "$new_archive"
symbols=$("$LLVM_NM" "$new_archive")
unresolved=$("$LLVM_NM" -u "$new_archive")
grep -q 'lfmm3d_t_cd_p_$' <<<"$symbols"
if grep -q 'lfmm3d_t_cd_p_c_$' <<<"$symbols"; then
  echo "fort2c archive unexpectedly contains differential-test symbols" >&2
  exit 1
fi
if grep -q '_c_$' <<<"$unresolved"; then
  echo "fort2c archive contains unresolved differential-test symbols" >&2
  exit 1
fi

# Publish only a fully generated and symbol-checked pair.
archive_next="$FMM3D_BUILD_DIR/.libfmm3d-wasm.a.next.$$"
stamp_next="$FMM3D_BUILD_DIR/.cache.stamp.next.$$"
cp "$new_archive" "$archive_next"
cp "$new_stamp" "$stamp_next"
mv -f "$archive_next" "$ARCHIVE"
mv -f "$stamp_next" "$STAMP"

echo "FMM3D_WASM_REBUILT archive=$ARCHIVE commit=$fmm_commit files=${#objects[@]}"
