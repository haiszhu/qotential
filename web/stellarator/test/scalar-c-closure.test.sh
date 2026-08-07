#!/usr/bin/env bash
# The scalar closure must lower cleanly both with and without the progress
# callback, and the callback symbol must appear only when it is enabled.
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "$0")/.." && pwd)
: "${LF:?set LF to the patched in-tree lfortran}"

out_plain=$(BIESOLVER_PROGRESS=0 LF="$LF" "$WEB_ROOT/scripts/probe-scalar-c-closure.sh")
printf '%s\n' "$out_plain"
grep -q '^SCALAR_C_CLOSURE_OK mode=plain$' <<<"$out_plain"

out_progress=$(BIESOLVER_PROGRESS=1 LF="$LF" "$WEB_ROOT/scripts/probe-scalar-c-closure.sh")
printf '%s\n' "$out_progress"
grep -q '^SCALAR_C_CLOSURE_OK mode=progress$' <<<"$out_progress"
