#!/usr/bin/env bash
# The 2-percent milestone helper must be deterministic and bounded.
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "$0")/.." && pwd)
TEST_DIR=$(mktemp -d "${TMPDIR:-/tmp}/stellarator-progress.XXXXXX")
trap 'rm -rf "$TEST_DIR"' EXIT

${FC:-gfortran} -cpp -O0 \
  "$WEB_ROOT/fortran/stellarator_progress_mod.f90" \
  "$WEB_ROOT/test/progress_milestone_test.f90" \
  -J "$TEST_DIR" -o "$TEST_DIR/progress-milestone"
"$TEST_DIR/progress-milestone"
