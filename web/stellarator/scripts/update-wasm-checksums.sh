#!/usr/bin/env bash
set -euo pipefail

PUBLIC_DIR=${PUBLIC_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../public/wasm" && pwd)}
files=(solver.js solver.wasm solver64.js solver64.wasm)
for file in "${files[@]}"; do
  test -s "$PUBLIC_DIR/$file" || {
    echo "missing or empty WASM artifact: $PUBLIC_DIR/$file" >&2
    exit 1
  }
done

next="$PUBLIC_DIR/.SHA256SUMS.next.$$"
trap 'rm -f "$next"' EXIT
(
  cd "$PUBLIC_DIR"
  shasum -a 256 "${files[@]}"
) > "$next"
test "$(wc -l < "$next" | tr -d ' ')" = 4
(
  cd "$PUBLIC_DIR"
  shasum -a 256 --check "$(basename "$next")" >/dev/null
)
mv -f "$next" "$PUBLIC_DIR/SHA256SUMS"
trap - EXIT
echo WASM_CHECKSUMS_UPDATED
