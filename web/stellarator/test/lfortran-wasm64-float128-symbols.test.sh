#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
LLVM_NM=${LLVM_NM:-llvm-nm}
OBJECT=${1:-$WEB_ROOT/build-wasm64/lfortran_intrinsics.o}
test -s "$OBJECT" || { echo "missing wasm64 runtime object: $OBJECT" >&2; exit 1; }

symbols=$("$LLVM_NM" "$OBJECT")
if grep -Eq '[[:space:]][TW][[:space:]]+__(add|sub|mul|div|neg)tf3$|[[:space:]][TW][[:space:]]+__(eq|ne|lt|le|gt|ge|unord)tf2$|[[:space:]][TW][[:space:]]+__(extend(df|sf)tf2|trunctf(df|sf)2|floatsi|floatdi|floatundi|fixtfsi|fixtfdi)$' \
    <<<"$symbols"; then
  echo "LFortran runtime exports compiler-rt binary128 symbols in wasm64" >&2
  grep -E '[[:space:]][TW][[:space:]]+__.*tf' <<<"$symbols" >&2 || true
  exit 1
fi
echo LFORTRAN_WASM64_FLOAT128_SYMBOLS_OK
