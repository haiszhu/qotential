#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "$0")/.." && pwd)
QOTENTIAL_ROOT=$(cd "$WEB_ROOT/../.." && pwd)
CORE="$WEB_ROOT/fortran/stellarator_grf_core_mod.f90"

test -s "$CORE"
grep -q '^module stellarator_grf_core_mod' "$CORE"
grep -q 'type, public :: stellarator_case_config' "$CORE"
grep -q 'type, public :: stellarator_case_result' "$CORE"
grep -q 'public :: stellarator_run_case, stellarator_result_clear' "$CORE"
grep -q 'use stellarator_grf_core_mod' "$QOTENTIAL_ROOT/test/stellarator_grf.f90"
grep -q 'call stellarator_run_case' "$QOTENTIAL_ROOT/test/stellarator_grf.f90"
! grep -q 'stellarator_run_cli' "$QOTENTIAL_ROOT/test/stellarator_grf.f90"

make -C "$QOTENTIAL_ROOT/test" sgrf >/tmp/stellarator-shared-core-build.log 2>&1
out=$(OMP_NUM_THREADS=1 "$QOTENTIAL_ROOT/build/stellarator_grf" 1 0 1 0)
printf '%s\n' "$out"
grep -q 'GRF max rel err =  *2.605E-04' <<<"$out"
