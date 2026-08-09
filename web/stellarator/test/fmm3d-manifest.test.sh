#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FMM3D_SRC=${FMM3D_SRC:-$HOME/git/FMM3D}
MANIFEST="$WEB_ROOT/fortran/fmm3d_wasm_sources.txt"

test -f "$MANIFEST"
count=0
while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  case "$relative" in
    src/Common/*.f|src/Laplace/*.f) ;;
    *) echo "unexpected FMM3D source: $relative" >&2; exit 1 ;;
  esac
  if [[ $relative == */prini.f ]]; then
    echo "diagnostic-only prini.f must not be in the WASM closure" >&2
    exit 1
  fi
  test -f "$FMM3D_SRC/$relative"
  count=$((count + 1))
done < "$MANIFEST"

test "$count" -eq 25
grep -qx 'src/Laplace/lfmm3dwrap.f' "$MANIFEST"
echo "FMM3D_MANIFEST_OK files=$count"
