#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
FORT2C=${FORT2C:-fort2c}
CLANG=${CLANG:-clang}
CLANG_CFLAGS_EXTRA=${CLANG_CFLAGS_EXTRA:-}
MANIFEST="$WEB_ROOT/fortran/fmm3d_wasm_sources.txt"
BUILD_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fmm3d-fort2c-native.XXXXXX")
trap 'rm -rf "$BUILD_DIR"' EXIT

mkdir -p "$BUILD_DIR/generated" "$BUILD_DIR/objects"
cp "$WEB_ROOT/native/fmm3d_c.h" "$BUILD_DIR/generated/fmm3d_c.h"

objects=()
while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  source="$FMM3D_SRC/$relative"
  name=$(basename "${relative%.f}")
  "$FORT2C" "$source" --runtime-header fmm3d_c.h \
    --guard-prefix BIESOLVER_FMM3D_ -o "$BUILD_DIR/generated" >/dev/null
  "$CLANG" -O2 -std=c99 -DFMM3D_DROP_IN -ffp-contract=off \
    -fno-strict-aliasing -fwrapv -Wno-unused-parameter \
    -Wno-unused-variable -Wno-parentheses -Wno-deprecated-non-prototype \
    $CLANG_CFLAGS_EXTRA \
    -I"$BUILD_DIR/generated" -c "$BUILD_DIR/generated/$name.c" \
    -o "$BUILD_DIR/objects/$name.o"
  objects+=("$BUILD_DIR/objects/$name.o")
done < "$MANIFEST"

"$CLANG" -O2 -std=c99 -ffp-contract=off \
  $CLANG_CFLAGS_EXTRA \
  "$SCRIPT_DIR/fmm3d_fort2c_native_test.c" "${objects[@]}" -lm \
  -o "$BUILD_DIR/fmm3d-fort2c-native-test"
"$BUILD_DIR/fmm3d-fort2c-native-test"
