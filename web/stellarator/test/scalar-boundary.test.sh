#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "$0")/.." && pwd)
QOTENTIAL_ROOT=$(cd "$WEB_ROOT/../.." && pwd)
QA_ROOT=$(cd "$QOTENTIAL_ROOT/external/QuatApproximation" && pwd)
LQ_ROOT=$(cd "$QOTENTIAL_ROOT/external/LineQuaaadrature" && pwd)
TMP=$(mktemp -d "${TMPDIR:-/tmp}/stellarator-scalar.XXXXXX")
trap 'rm -rf "$TMP"' EXIT

FC=${FC:-gfortran-16}
FLAGS=(-O0 -cpp -DBIESOLVER_WASM_SCALAR_ONLY -DBIESOLVER_R64_ONLY \
  -DBIESOLVER_STELLARATOR_BUILD -fdefault-integer-8 -J"$TMP" -I"$TMP" \
  -I"$QOTENTIAL_ROOT/build" -I"$QA_ROOT/build" -I"$LQ_ROOT/build")
"$FC" "${FLAGS[@]}" -c "$LQ_ROOT/src/kernel_mod.f90" -o "$TMP/kernel.o"
"$FC" "${FLAGS[@]}" -c "$LQ_ROOT/src/solidangle_mod.f90" -o "$TMP/solidangle.o"
"$FC" "${FLAGS[@]}" -c "$QOTENTIAL_ROOT/src/patch_refine_mod.f90" -o "$TMP/patch.o"
"$FC" "${FLAGS[@]}" -c "$WEB_ROOT/fortran/w7x_modes_mod.f90" -o "$TMP/w7x.o"
"$FC" "${FLAGS[@]}" -c "$WEB_ROOT/fortran/stellarator_grf_core_mod.f90" \
  -o "$TMP/core.o"
test -s "$TMP/kernel.o" -a -s "$TMP/solidangle.o" -a -s "$TMP/patch.o" -a -s "$TMP/core.o"

bad=$(nm -u "$TMP"/*.o | grep -Ei 'lfmm3d|omp_get_|__kmpc|GOMP_|lq_csimd|csimd|cpu_time|system_clock|sys_clock' || true)
if [[ -n $bad ]]; then
  printf '%s\n' "$bad" >&2
  echo 'scalar closure still references FMM/OpenMP/SIMD/native clocks' >&2
  exit 1
fi

"$FC" "${FLAGS[@]}" -DBIESOLVER_WASM_FMM3D \
  -c "$WEB_ROOT/fortran/stellarator_grf_core_mod.f90" -o "$TMP/core-fmm.o"
fmm_symbols=$(nm -u "$TMP/core-fmm.o")
adapter_count=$(grep -Ec 'biesolver_lfmm3d_t_cd_p_rowmajor' <<<"$fmm_symbols" || true)
raw_count=$(grep -Ec 'lfmm3d_t_cd_p_$' <<<"$fmm_symbols" || true)
if [[ $adapter_count != 1 || $raw_count != 0 ]]; then
  printf '%s\n' "$fmm_symbols" >&2
  echo "scalar+FMM closure must expose only the row-major adapter seam" >&2
  exit 1
fi
