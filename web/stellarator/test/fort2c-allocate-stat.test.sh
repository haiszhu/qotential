#!/usr/bin/env bash
set -euo pipefail

FORT2C=${FORT2C:-fort2c}
WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/fort2c-allocate-stat.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT

"$FORT2C" "$WEB_ROOT/test/fixtures/fort2c_alloc_probe.f" \
  --basename alloc_probe --runtime-header fmm3d_c.h \
  --guard-prefix BIESOLVER_FMM3D_ -o "$TMP_DIR" >/dev/null

grep -Eq 'ier[[:space:]]*=[[:space:]]*\([^;]*a[[:space:]]*==[[:space:]]*NULL[^;]*\)[[:space:]]*\?[[:space:]]*1[[:space:]]*:[[:space:]]*0;' \
  "$TMP_DIR/alloc_probe.c"
echo "FORT2C_ALLOCATE_STAT_OK"
