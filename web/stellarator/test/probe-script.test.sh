#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
: "${LF:?set LF to the patched in-tree lfortran}"

out=$(mktemp)
trap 'rm -f "$out"' EXIT

LF="$LF" "$ROOT/scripts/probe-lfortran-c-backend.sh" | tee "$out"
grep -Eq '^lap3d_close[[:space:]]+0[[:space:]]+[1-9][0-9]{3,}[[:space:]]+OK$' "$out"
grep -Eq '^stellarator_wasm_api[[:space:]]+0[[:space:]]+[1-9][0-9]{3,}[[:space:]]+OK$' "$out"
grep -q '^C_AND_WASM_OBJECT_GATE_OK$' "$out"
