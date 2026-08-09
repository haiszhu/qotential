#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
QOTENTIAL_DIR=$(cd "$SCRIPT_DIR/../../.." && pwd)
QA_DIR="$QOTENTIAL_DIR/external/QuatApproximation"
LQ_DIR="$QOTENTIAL_DIR/external/LineQuaaadrature"
LF=${LF:-lfortran}
EMCC=${EMCC:-emcc}
LLVM_NM=${LLVM_NM:-llvm-nm}

LF_BIN=$(command -v "$LF" 2>/dev/null || printf '%s\n' "$LF")
if [[ -z ${LF_SRC:-} ]]; then
  LF_SRC=$(cd "$(dirname "$LF_BIN")/../../.." && pwd)
fi
LF_RUNTIME_DIR="$LF_SRC/src/libasr/runtime"
test -f "$LF_RUNTIME_DIR/lfortran_intrinsics.h"

PROBE_DIR=$(mktemp -d "${TMPDIR:-/tmp}/stellarator-scalar-c.XXXXXX")
trap 'rm -rf "$PROBE_DIR"' EXIT
MOD_DIR="$PROBE_DIR/modules"
OBJ_DIR="$PROBE_DIR/c-objects"
mkdir -p "$MOD_DIR" "$OBJ_DIR"

CPP_FLAGS=(
  --cpp --no-style-suggestions --legacy-array-sections
  --no-array-bounds-checking
  -D BIESOLVER_R64_ONLY
  -D BIESOLVER_STELLARATOR_BUILD
  -D BIESOLVER_WASM_SCALAR_ONLY
  -D BIESOLVER_C_BACKEND_ROW_MAJOR
)
# Progress instrumentation is opt-in so the same closure can be probed with the
# callback present and absent.
BIESOLVER_PROGRESS=${BIESOLVER_PROGRESS:-0}
if [[ $BIESOLVER_PROGRESS == 1 ]]; then
  CPP_FLAGS+=(-D BIESOLVER_WASM_PROGRESS)
fi
BIESOLVER_FMM3D=${BIESOLVER_FMM3D:-0}
if [[ $BIESOLVER_FMM3D == 1 ]]; then
  CPP_FLAGS+=(-D BIESOLVER_WASM_FMM3D)
fi
NAMES=(
  quatapproximation harmonic koorn_geom qkernel omega tensor_geom
  linequaaadrature adaptive kernel solidangle lap3d
  patch_refine lap3d_close w7x_modes stellarator_progress stellarator_grf_core
)
SOURCES=(
  "$QA_DIR/src/quatapproximation_mod.f90"
  "$QA_DIR/src/harmonic_mod.f90"
  "$QA_DIR/src/koorn_geom_mod.f90"
  "$QA_DIR/src/qkernel_mod.f90"
  "$QA_DIR/src/omega_mod.f90"
  "$QA_DIR/src/tensor_geom_mod.f90"
  "$LQ_DIR/src/linequaaadrature_mod.f90"
  "$LQ_DIR/src/adaptive_mod.f90"
  "$LQ_DIR/src/kernel_mod.f90"
  "$LQ_DIR/src/solidangle_mod.f90"
  "$LQ_DIR/src/lap3d_mod.f90"
  "$QOTENTIAL_DIR/src/patch_refine_mod.f90"
  "$QOTENTIAL_DIR/src/lap3d_close_mod.f90"
  "$SCRIPT_DIR/../fortran/w7x_modes_mod.f90"
  "$SCRIPT_DIR/../fortran/stellarator_progress_mod.f90"
  "$SCRIPT_DIR/../fortran/stellarator_grf_core_mod.f90"
)

for i in "${!NAMES[@]}"; do
  name=${NAMES[$i]}
  source=${SOURCES[$i]}
  echo "SETUP $name"
  "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
    -c "$source" -o "$PROBE_DIR/mod-$name.o"
  test -s "$PROBE_DIR/mod-$name.o"
done

for i in "${!NAMES[@]}"; do
  name=${NAMES[$i]}
  source=${SOURCES[$i]}
  echo "C $name"
  c_file="$PROBE_DIR/$name.c"
  object="$OBJ_DIR/$name.o"
  "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
    --show-c "$source" >"$c_file"
  lines=$(wc -l <"$c_file" | tr -d ' ')
  if (( lines < 20 )); then
    echo "$name: generated C is only $lines lines" >&2
    exit 1
  fi
  "$EMCC" -std=c11 -ffp-contract=off -I"$LF_RUNTIME_DIR" -c "$c_file" -o "$object"
  test -s "$object"
done

bad=''
for object in "$OBJ_DIR"/*.o; do
  matches=$("$LLVM_NM" -u "$object" | \
    grep -Ei 'lq_csimd|csimd|cpu_time|system_clock|sys_clock|omp_get_|__kmpc|GOMP_' || true)
  if [[ -n $matches ]]; then
    bad+="$(basename "$object")"$'\n'"$matches"$'\n'
  fi
done
if [[ -n $bad ]]; then
  printf '%s\n' "$bad" >&2
  exit 1
fi

core_undef=$($LLVM_NM -u "$OBJ_DIR/stellarator_grf_core.o")
adapter_count=$(grep -Ec '(^|[[:space:]])biesolver_lfmm3d_t_cd_p_rowmajor$' \
  <<<"$core_undef" || true)
raw_fmm_count=$(grep -Ec '(^|[[:space:]])lfmm3d_t_cd_p_$' <<<"$core_undef" || true)
if [[ $BIESOLVER_FMM3D == 1 ]]; then
  if [[ $adapter_count != 1 || $raw_fmm_count != 0 ]]; then
    echo "stellarator_grf_core.o: scalar+FMM must expose only the row-major adapter" >&2
    exit 1
  fi
elif [[ $adapter_count != 0 || $raw_fmm_count != 0 ]]; then
  echo "stellarator_grf_core.o: scalar Direct closure unexpectedly references FMM" >&2
  exit 1
fi
preflight_count=$(grep -Ec '(^|[[:space:]])biesolver_wasm_memory_preflight$' <<<"$core_undef" || true)
if [[ $preflight_count != 1 ]]; then
  echo "stellarator_grf_core.o: expected exactly one undefined biesolver_wasm_memory_preflight, found $preflight_count" >&2
  exit 1
fi
for symbol in biesolver_charts_rowmajor biesolver_geo_charts_rowmajor \
              biesolver_uv2x_charts_rowmajor; do
  if ! grep -q "$symbol" <<<"$core_undef"; then
    echo "stellarator_grf_core.o: missing generated call to $symbol" >&2
    exit 1
  fi
done
if grep -Eq '(^|[[:space:]])stellarator_(charts|geo_charts|uv2x_charts)$' \
     <<<"$core_undef"; then
  echo "stellarator_grf_core.o: bypasses the row-major geometry adapter" >&2
  exit 1
fi
progress_count=$(grep -Ec '(^|[[:space:]])biesolver_progress$' <<<"$core_undef" || true)
if [[ $BIESOLVER_PROGRESS == 1 ]]; then
  if [[ $progress_count != 1 ]]; then
    echo "stellarator_grf_core.o: expected exactly one undefined biesolver_progress, found $progress_count" >&2
    exit 1
  fi
  echo "SCALAR_C_CLOSURE_OK mode=progress fmm=$BIESOLVER_FMM3D"
else
  if [[ $progress_count != 0 ]]; then
    echo "stellarator_grf_core.o: biesolver_progress must be absent without the macro, found $progress_count" >&2
    exit 1
  fi
  echo "SCALAR_C_CLOSURE_OK mode=plain fmm=$BIESOLVER_FMM3D"
fi
