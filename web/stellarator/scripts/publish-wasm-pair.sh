#!/usr/bin/env bash
set -euo pipefail

if [[ $# != 4 ]]; then
  echo "usage: publish-wasm-pair.sh <staged.js> <staged.wasm> <public-dir> <basename>" >&2
  exit 2
fi
STAGED_JS=$1
STAGED_WASM=$2
PUBLIC_DIR=$3
BASENAME=$4
if [[ $BASENAME != solver && $BASENAME != solver64 ]]; then
  echo "unsupported solver basename: $BASENAME" >&2
  exit 2
fi
test -s "$STAGED_JS" || { echo "missing or empty staged JavaScript: $STAGED_JS" >&2; exit 1; }
test -s "$STAGED_WASM" || { echo "missing or empty staged WASM: $STAGED_WASM" >&2; exit 1; }
mkdir -p "$PUBLIC_DIR"

if [[ ${WASM_PUBLISH_FAIL_BEFORE_RENAME:-0} == 1 ]]; then
  echo "injected failure before WASM publication" >&2
  exit 99
fi

# Publish the binary first and the loader last. Both staged members have already
# passed validation; a failed build never reaches this point.
mv -f "$STAGED_WASM" "$PUBLIC_DIR/$BASENAME.wasm"
mv -f "$STAGED_JS" "$PUBLIC_DIR/$BASENAME.js"
echo "WASM_PAIR_PUBLISHED basename=$BASENAME"
