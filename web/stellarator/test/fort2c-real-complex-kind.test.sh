#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FORT2C=${FORT2C:-fort2c}
FORT2C_SRC=${FORT2C_SRC:-$HOME/fort2c}
CLANG=${CLANG:-clang}
GFORTRAN=${GFORTRAN:-gfortran}
BASE_COMMIT=7f7fda827260df4f7ff1bfcaf1c420e3e809ac7b
FORT2C_CHECK_BASELINE=${FORT2C_CHECK_BASELINE:-1}
FIXTURE="$WEB_ROOT/test/fixtures/fort2c_real_complex_kind.f"
RUNTIME_HEADER="$WEB_ROOT/native/fmm3d_c.h"
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fort2c-real-complex-kind.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

for tool in "$FORT2C" "$CLANG" "$GFORTRAN" git python3 tar; do
  command -v "$tool" >/dev/null 2>&1 || {
    echo "required tool not found: $tool" >&2
    exit 1
  }
done
if [[ $FORT2C_CHECK_BASELINE != 0 && $FORT2C_CHECK_BASELINE != 1 ]]; then
  echo "FORT2C_CHECK_BASELINE must be 0 or 1" >&2
  exit 1
fi

if [[ $FORT2C_CHECK_BASELINE == 1 ]]; then
  test -d "$FORT2C_SRC/.git" || {
    echo "FORT2C_SRC is not a git checkout: $FORT2C_SRC" >&2
    exit 1
  }
  FORT2C_PATH=$(command -v "$FORT2C")
  FORT2C_REAL=$(python3 - "$FORT2C_PATH" <<'PY'
import os
import sys
print(os.path.realpath(sys.argv[1]))
PY
  )
  FORT2C_PYTHON=$(sed -n '1s/^#!//p' "$FORT2C_REAL")
  test -x "$FORT2C_PYTHON" || {
    echo "cannot resolve fort2c Python interpreter from $FORT2C_REAL" >&2
    exit 1
  }
fi

mkdir -p "$TMP_DIR/fixed"
cp "$RUNTIME_HEADER" "$TMP_DIR/fixed/fmm3d_c.h"

"$FORT2C" "$FIXTURE" --runtime-header fmm3d_c.h \
  -o "$TMP_DIR/fixed" >/dev/null

FIXED_C="$TMP_DIR/fixed/fort2c_real_complex_kind.c"
if grep -Fq '(float)creal(' "$FIXED_C"; then
  echo "corrected fort2c still narrows REAL(COMPLEX*16) to float" >&2
  exit 1
fi

"$CLANG" -O0 -std=c99 "$FIXED_C" -lm -o "$TMP_DIR/fixed.exe"
"$GFORTRAN" -O0 "$FIXTURE" -o "$TMP_DIR/reference.exe"

"$TMP_DIR/fixed.exe" > "$TMP_DIR/fixed.out"
"$TMP_DIR/reference.exe" > "$TMP_DIR/reference.out"
cmp "$TMP_DIR/reference.out" "$TMP_DIR/fixed.out"

if [[ $FORT2C_CHECK_BASELINE == 1 ]]; then
  mkdir -p "$TMP_DIR/base-src" "$TMP_DIR/base"
  git -C "$FORT2C_SRC" archive "$BASE_COMMIT" | tar -x -C "$TMP_DIR/base-src"
  cp "$RUNTIME_HEADER" "$TMP_DIR/base/fmm3d_c.h"
  PYTHONPATH="$TMP_DIR/base-src" "$FORT2C_PYTHON" -m fort2c.cli \
    "$FIXTURE" --runtime-header fmm3d_c.h -o "$TMP_DIR/base" >/dev/null
  BASE_C="$TMP_DIR/base/fort2c_real_complex_kind.c"
  grep -Fq '(float)creal(' "$BASE_C"
  "$CLANG" -O0 -std=c99 "$BASE_C" -lm -o "$TMP_DIR/base.exe"
  "$TMP_DIR/base.exe" > "$TMP_DIR/base.out"
  if cmp -s "$TMP_DIR/base.out" "$TMP_DIR/reference.out"; then
    echo "pre-fix fort2c unexpectedly matches gfortran" >&2
    exit 1
  fi
  printf 'FORT2C_REAL_COMPLEX_KIND_OK base=%s fixed=%s\n' \
    "$(tr -d ' ' < "$TMP_DIR/base.out")" \
    "$(tr -d ' ' < "$TMP_DIR/fixed.out")"
else
  printf 'FORT2C_REAL_COMPLEX_KIND_OK fixed=%s\n' \
    "$(tr -d ' ' < "$TMP_DIR/fixed.out")"
fi
