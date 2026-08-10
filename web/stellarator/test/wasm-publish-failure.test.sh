#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
PUBLISHER="$WEB_ROOT/scripts/publish-wasm-pair.sh"
TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/wasm-publish.XXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT
PUBLIC_DIR="$TMP_DIR/public"
STAGE_DIR="$TMP_DIR/stage"
mkdir -p "$PUBLIC_DIR" "$STAGE_DIR"

printf 'old-js\n' > "$PUBLIC_DIR/solver64.js"
printf 'old-wasm\n' > "$PUBLIC_DIR/solver64.wasm"
printf 'old-checksum\n' > "$PUBLIC_DIR/SHA256SUMS"
cp "$PUBLIC_DIR/solver64.js" "$TMP_DIR/expected.js"
cp "$PUBLIC_DIR/solver64.wasm" "$TMP_DIR/expected.wasm"
cp "$PUBLIC_DIR/SHA256SUMS" "$TMP_DIR/expected.sha"

printf 'new-js\n' > "$STAGE_DIR/solver64.js"
: > "$STAGE_DIR/solver64.wasm"
if "$PUBLISHER" "$STAGE_DIR/solver64.js" "$STAGE_DIR/solver64.wasm" \
    "$PUBLIC_DIR" solver64; then
  echo 'publisher accepted an empty WASM member' >&2
  exit 1
fi
cmp "$TMP_DIR/expected.js" "$PUBLIC_DIR/solver64.js"
cmp "$TMP_DIR/expected.wasm" "$PUBLIC_DIR/solver64.wasm"
cmp "$TMP_DIR/expected.sha" "$PUBLIC_DIR/SHA256SUMS"

printf 'new-wasm\n' > "$STAGE_DIR/solver64.wasm"
if WASM_PUBLISH_FAIL_BEFORE_RENAME=1 \
    "$PUBLISHER" "$STAGE_DIR/solver64.js" "$STAGE_DIR/solver64.wasm" \
    "$PUBLIC_DIR" solver64; then
  echo 'injected pre-publication failure unexpectedly succeeded' >&2
  exit 1
fi
cmp "$TMP_DIR/expected.js" "$PUBLIC_DIR/solver64.js"
cmp "$TMP_DIR/expected.wasm" "$PUBLIC_DIR/solver64.wasm"
cmp "$TMP_DIR/expected.sha" "$PUBLIC_DIR/SHA256SUMS"

"$PUBLISHER" "$STAGE_DIR/solver64.js" "$STAGE_DIR/solver64.wasm" \
  "$PUBLIC_DIR" solver64
grep -q '^new-js$' "$PUBLIC_DIR/solver64.js"
grep -q '^new-wasm$' "$PUBLIC_DIR/solver64.wasm"
cmp "$TMP_DIR/expected.sha" "$PUBLIC_DIR/SHA256SUMS"
echo WASM_PUBLISH_FAILURE_PRESERVATION_OK
