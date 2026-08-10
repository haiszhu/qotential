#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
WEB_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
QOTENTIAL_DIR=$(cd "$WEB_ROOT/../.." && pwd)
LF=${LF:-lfortran}
EMCC=${EMCC:-emcc}
LLVM_NM=${LLVM_NM:-llvm-nm}
EMAR=${EMAR:-emar}
WASM_OPT=${WASM_OPT:--O2}
WASM_DEBUG=${WASM_DEBUG:-}
WASM_ALLOC_TRACE=${WASM_ALLOC_TRACE:-0}
LF_ARRAY_BOUNDS_CHECK=${LF_ARRAY_BOUNDS_CHECK:-0}
WASM_SANITIZE=${WASM_SANITIZE:-}
WASM_STACK_SIZE=${WASM_STACK_SIZE:-67108864}
WASM_PROGRESS=${WASM_PROGRESS:-1}
WASM_MEMORY64=${WASM_MEMORY64:-0}

sha256_file() {
  shasum -a 256 "$1" | awk '{print $1}'
}
if [[ $WASM_PROGRESS != 0 && $WASM_PROGRESS != 1 ]]; then
  echo "WASM_PROGRESS must be 0 or 1 (got '$WASM_PROGRESS')" >&2
  exit 1
fi
if [[ $WASM_MEMORY64 != 0 && $WASM_MEMORY64 != 1 ]]; then
  echo "WASM_MEMORY64 must be 0 or 1 (got '$WASM_MEMORY64')" >&2
  exit 1
fi
LF_BIN=$(command -v "$LF" 2>/dev/null || printf '%s\n' "$LF")
if [[ ! -x $LF_BIN ]]; then
  echo "LF=$LF is not an executable." >&2
  echo "Build the patched compiler first; see the macOS toolchain section of" >&2
  echo "web/stellarator/README.md.  LF must be <lfortran-source>/build/src/bin/lfortran" >&2
  echo "-- the in-tree binary, because the source tree is recovered from its path." >&2
  exit 1
fi
if [[ -z ${LF_SRC:-} ]]; then
  LF_SRC=$(cd "$(dirname "$LF_BIN")/../../.." && pwd)
fi
RUNTIME_DIR="$LF_SRC/src/libasr/runtime"
LFORTRAN_WASM64_PATCH="$WEB_ROOT/patches/lfortran-wasm64-float128-symbols.patch"
if [[ ! -f $RUNTIME_DIR/lfortran_intrinsics.h || ! -f $RUNTIME_DIR/lfortran_intrinsics.c ]]; then
  echo "no LFortran runtime headers under $RUNTIME_DIR" >&2
  echo "LF_SRC was derived as $LF_SRC from LF=$LF_BIN; it must be the LFortran" >&2
  echo "source tree.  Point LF at <source>/build/src/bin/lfortran, or set LF_SRC." >&2
  exit 1
fi
test -f "$LFORTRAN_WASM64_PATCH" || {
  echo "missing LFortran wasm64 runtime patch: $LFORTRAN_WASM64_PATCH" >&2
  exit 1
}

if [[ $WASM_MEMORY64 == 1 ]]; then
  architecture=wasm64
  BUILD_DIR="$WEB_ROOT/build-wasm64"
  FMM3D_BUILD_DIR="$WEB_ROOT/build-fmm3d-wasm64"
  OUTPUT_BASENAME=solver64
else
  architecture=wasm32
  BUILD_DIR="$WEB_ROOT/build-wasm32"
  FMM3D_BUILD_DIR="$WEB_ROOT/build-fmm3d-wasm32"
  OUTPUT_BASENAME=solver
fi
MOD_DIR="$BUILD_DIR/modules"
PUBLIC_DIR=${PUBLIC_DIR:-$WEB_ROOT/public/wasm}
rm -rf "$BUILD_DIR"
mkdir -p "$BUILD_DIR" "$MOD_DIR" "$PUBLIC_DIR"

EMCC_FLAGS=("$WASM_OPT")
if [[ -n $WASM_DEBUG ]]; then
  EMCC_FLAGS+=("$WASM_DEBUG")
fi
if [[ $WASM_MEMORY64 == 1 ]]; then
  EMCC_FLAGS+=(-m64)
fi

CPP_FLAGS=(
  --cpp --no-style-suggestions --legacy-array-sections
  -D BIESOLVER_R64_ONLY
  -D BIESOLVER_STELLARATOR_BUILD
  -D BIESOLVER_WASM_SCALAR_ONLY
  -D BIESOLVER_WASM_FMM3D
  -D BIESOLVER_C_BACKEND_ROW_MAJOR
)
if [[ $LF_ARRAY_BOUNDS_CHECK == 0 ]]; then
  CPP_FLAGS+=(--no-array-bounds-checking)
fi
# Progress instrumentation is opt-in and lives entirely behind this macro.  With
# WASM_PROGRESS=0 neither the macro nor the EM_JS bridge object is present, so
# the generated C and solver.wasm reproduce the pre-progress baseline byte-for-byte.
if [[ $WASM_PROGRESS == 1 ]]; then
  CPP_FLAGS+=(-D BIESOLVER_WASM_PROGRESS)
fi
SANITIZER_FLAGS=(-DBIESOLVER_WASM_BUILD=1)
PATH_MAP="-ffile-prefix-map=$WEB_ROOT=."
if [[ -n $WASM_SANITIZE ]]; then
  SANITIZER_FLAGS+=("-fsanitize=$WASM_SANITIZE")
fi

while IFS= read -r relative; do
  [[ -z $relative || $relative == \#* ]] && continue
  source="$QOTENTIAL_DIR/$relative"
  name=$(basename "${relative%.f90}")
  "$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" \
    -c "$source" -o "$BUILD_DIR/setup-$name.o"
  test -s "$BUILD_DIR/setup-$name.o"
done < "$WEB_ROOT/fortran/wasm_sources.txt"

SOLVER_C="$BUILD_DIR/stellarator_solver.c"
"$LF" "${CPP_FLAGS[@]}" -J "$MOD_DIR" -I "$MOD_DIR" --show-c \
  "$WEB_ROOT/fortran/stellarator_grf_core_mod.f90" > "$SOLVER_C"
python3 "$WEB_ROOT/scripts/rewrite-lfortran-c.py" "$SOLVER_C"
lines=$(wc -l < "$SOLVER_C" | tr -d ' ')
if (( lines < 1000 )); then
  echo "generated solver C is only $lines lines" >&2
  exit 1
fi

# Not under /tmp: on macOS that is a symlink to /private/tmp, and Emscripten
# builds its cached system libraries with source paths computed relative to the
# cache directory.  The symlink throws that arithmetic off by a level and the
# build dies with "clang: error: no such file or directory:
# '../../../../opt/homebrew/Cellar/emscripten/.../gl.c'".  Any real path works;
# `em-config CACHE` reports Emscripten's own default if you prefer that.
export EM_CACHE=${EM_CACHE:-${XDG_CACHE_HOME:-$HOME/.cache}/biesolver-emscripten}
FMM3D_CFLAGS_EXTRA=
if [[ -n $WASM_SANITIZE ]]; then
  FMM3D_CFLAGS_EXTRA="-fsanitize=$WASM_SANITIZE"
fi
FMM3D_SRC="${FMM3D_SRC:-$HOME/git/FMM3D}" FORT2C="${FORT2C:-fort2c}" \
  FORT2C_SRC="${FORT2C_SRC:-$HOME/fort2c}" \
  EMCC="$EMCC" EMAR="$EMAR" \
  WASM_MEMORY64="$WASM_MEMORY64" FMM3D_BUILD_DIR="$FMM3D_BUILD_DIR" \
  FMM3D_CFLAGS_EXTRA="$FMM3D_CFLAGS_EXTRA" \
  "$WEB_ROOT/scripts/build-fmm3d-wasm.sh"
SOLVER_CC_FLAGS=(
  -DBIESOLVER_GENERATED_SOLVER=1
  -include "$WEB_ROOT/native/wasm_lfortran_alloc.h"
)
if [[ $WASM_ALLOC_TRACE == 1 ]]; then
  SOLVER_CC_FLAGS+=(
    -Dmalloc=biesolver_debug_malloc
    -Dfree=biesolver_debug_free
  )
fi
"$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 -ffp-contract=off -I"$RUNTIME_DIR" \
  "${SANITIZER_FLAGS[@]}" \
  "${SOLVER_CC_FLAGS[@]}" \
  -c "$SOLVER_C" -o "$BUILD_DIR/stellarator_solver.o"
test -s "$BUILD_DIR/stellarator_solver.o"

undef=$($LLVM_NM -u "$BUILD_DIR/stellarator_solver.o")
adapter_count=$(grep -Ec '(^|[[:space:]])biesolver_lfmm3d_t_cd_p_rowmajor$' \
  <<<"$undef" || true)
raw_fmm_count=$(grep -Ec '(^|[[:space:]])lfmm3d_t_cd_p_$' <<<"$undef" || true)
if [[ $adapter_count != 1 || $raw_fmm_count != 0 ]]; then
  echo "generated solver must reference only the FMM3D row-major adapter" >&2
  exit 1
fi
if grep -Ei 'lq_csimd|csimd|omp_get_|__kmpc|GOMP_|cpu_time|system_clock' \
     <<<"$undef"; then
  echo "forbidden scalar-build dependency remains" >&2
  exit 1
fi

"$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 -ffp-contract=off \
  "${SANITIZER_FLAGS[@]}" -I"$WEB_ROOT/native" \
  -c "$WEB_ROOT/native/fmm3d_layout_adapter.c" \
  -o "$BUILD_DIR/fmm3d_layout_adapter.o"
test -s "$BUILD_DIR/fmm3d_layout_adapter.o"
adapter_undef=$($LLVM_NM -u "$BUILD_DIR/fmm3d_layout_adapter.o")
adapter_raw_count=$(grep -Ec '(^|[[:space:]])lfmm3d_t_cd_p_$' \
  <<<"$adapter_undef" || true)
if [[ $adapter_raw_count != 1 ]]; then
  echo "FMM3D layout adapter must reference the canonical raw FMM symbol" >&2
  exit 1
fi

RUNTIME_BUILD_DIR="$BUILD_DIR/lfortran-runtime"
mkdir -p "$RUNTIME_BUILD_DIR"
cp "$RUNTIME_DIR/lfortran_intrinsics.c" "$RUNTIME_BUILD_DIR/"
cp "$RUNTIME_DIR/lfortran_float128.c" "$RUNTIME_BUILD_DIR/"
cp "$RUNTIME_DIR/lfortran_float128_llvm.h" "$RUNTIME_BUILD_DIR/"
patch -s -d "$RUNTIME_BUILD_DIR" -p1 < "$LFORTRAN_WASM64_PATCH"

"$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 -D_GNU_SOURCE -I"$LF_SRC/src" -I"$RUNTIME_DIR" \
  "${SANITIZER_FLAGS[@]}" \
  -D_lfortran_get_default_allocator=_lfortran_runtime_get_default_allocator \
  -D_lfortran_malloc_alloc=_lfortran_runtime_malloc_alloc \
  -D_lfortran_free_alloc=_lfortran_runtime_free_alloc \
  -c "$RUNTIME_BUILD_DIR/lfortran_intrinsics.c" -o "$BUILD_DIR/lfortran_intrinsics.o"
"$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 "${SANITIZER_FLAGS[@]}" \
  -c "$WEB_ROOT/native/wasm_lfortran_alloc.c" \
  -o "$BUILD_DIR/wasm_lfortran_alloc.o"

# Under WASM_PROGRESS=1 the EM_JS bridge resolves the Fortran core's
# `biesolver_progress` import at the final link.  The callback stays an internal
# import; it is never added to EXPORTED_FUNCTIONS.
PROGRESS_LINK_OBJS=()
if [[ $WASM_PROGRESS == 1 ]]; then
  "$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 "${SANITIZER_FLAGS[@]}" \
    -c "$WEB_ROOT/native/wasm_progress_bridge.c" \
    -o "$BUILD_DIR/wasm_progress_bridge.o"
  test -s "$BUILD_DIR/wasm_progress_bridge.o"
  PROGRESS_LINK_OBJS=("$BUILD_DIR/wasm_progress_bridge.o")
fi

EXPORTED="['_solver_simplex_precomp','_solver_run','_solver_result_nsrc','_solver_result_nrender','_solver_result_ntriangles','_solver_result_grf_error','_solver_copy_sx','_solver_copy_snx','_solver_copy_sw','_solver_copy_ub','_solver_copy_ubn','_solver_copy_u','_solver_copy_render_xyz','_solver_copy_render_log_error','_solver_copy_render_triangles','_solver_last_error','_solver_clear','_malloc','_free']"

MEMORY_LINK_FLAGS=(-sALLOW_MEMORY_GROWTH=1)
if [[ $WASM_MEMORY64 == 1 ]]; then
  MEMORY_LINK_FLAGS+=(-sINITIAL_MEMORY=134217728 -sMAXIMUM_MEMORY=17179869184)
fi
PUBLISH_STAGE="$BUILD_DIR/publish"
mkdir -p "$PUBLISH_STAGE"

"$EMCC" "${EMCC_FLAGS[@]}" "$PATH_MAP" -std=c11 -D_GNU_SOURCE -ffp-contract=off \
  "${SANITIZER_FLAGS[@]}" \
  -Wno-unused-variable -I"$LF_SRC/src" -I"$RUNTIME_DIR" -I"$WEB_ROOT/native" \
  "$BUILD_DIR/stellarator_solver.o" \
  "$BUILD_DIR/lfortran_intrinsics.o" \
  "$BUILD_DIR/wasm_lfortran_alloc.o" \
  "$BUILD_DIR/fmm3d_layout_adapter.o" \
  ${PROGRESS_LINK_OBJS[@]+"${PROGRESS_LINK_OBJS[@]}"} \
  "$QOTENTIAL_DIR/utils/stellarator_geo_mex.c" \
  "$WEB_ROOT/native/geometry_layout_adapter.c" \
  "$WEB_ROOT/native/wasm_blas_shim.c" \
  "$WEB_ROOT/native/wasm_memory_preflight.c" \
  "$WEB_ROOT/native/wasm_api_adapter.c" \
  "$FMM3D_BUILD_DIR/libfmm3d-wasm.a" \
  -sMODULARIZE=1 -sEXPORT_ES6=1 -sEXPORT_NAME=createSolver \
  -sWASM_BIGINT=1 "${MEMORY_LINK_FLAGS[@]}" -sENVIRONMENT=web,worker,node \
  -sSTACK_SIZE="$WASM_STACK_SIZE" \
  -sINCOMING_MODULE_JS_API=wasmBinary,locateFile \
  -sEXPORTED_RUNTIME_METHODS=HEAPU8 \
  -sEXPORTED_FUNCTIONS="$EXPORTED" \
  -o "$PUBLISH_STAGE/$OUTPUT_BASENAME.js" -lm

test -s "$PUBLISH_STAGE/$OUTPUT_BASENAME.js"
test -s "$PUBLISH_STAGE/$OUTPUT_BASENAME.wasm"
node --input-type=module -e \
  'const m=await import(process.argv[1]); await m.default();' \
  "file://$PUBLISH_STAGE/$OUTPUT_BASENAME.js"
"$WEB_ROOT/scripts/publish-wasm-pair.sh" \
  "$PUBLISH_STAGE/$OUTPUT_BASENAME.js" \
  "$PUBLISH_STAGE/$OUTPUT_BASENAME.wasm" "$PUBLIC_DIR" "$OUTPUT_BASENAME"
test -s "$PUBLIC_DIR/$OUTPUT_BASENAME.js"
test -s "$PUBLIC_DIR/$OUTPUT_BASENAME.wasm"

FMM3D_STAMP="$FMM3D_BUILD_DIR/cache.stamp"
test -s "$FMM3D_STAMP"
lfortran_commit=$(git -C "$LF_SRC" rev-parse HEAD)
lfortran_describe=$(git -C "$LF_SRC" describe --tags --always --dirty)
lfortran_diff_sha=$(git -C "$LF_SRC" diff --binary HEAD | shasum -a 256 | awk '{print $1}')
provenance_next="$BUILD_DIR/.provenance.stamp.next.$$"
{
  printf 'architecture=%s\n' "$architecture"
  printf 'lfortran_commit=%s\n' "$lfortran_commit"
  printf 'lfortran_describe=%s\n' "$lfortran_describe"
  printf 'lfortran_diff_sha256=%s\n' "$lfortran_diff_sha"
  printf 'lfortran_binary_sha256=%s\n' "$(sha256_file "$LF_BIN")"
  printf 'lfortran_runtime_patch_sha256=%s\n' "$(sha256_file "$LFORTRAN_WASM64_PATCH")"
  printf 'generated_solver_c_sha256=%s\n' "$(sha256_file "$SOLVER_C")"
  printf 'emcc_version=%s\n' "$("$EMCC" --version | sed -n '1p')"
  printf 'emcc_flags='
  separator=
  for flag in "${EMCC_FLAGS[@]}" "${MEMORY_LINK_FLAGS[@]}"; do
    printf '%s%q' "$separator" "$flag"
    separator=' '
  done
  printf '\n'
  printf 'artifact_js_sha256=%s\n' "$(sha256_file "$PUBLIC_DIR/$OUTPUT_BASENAME.js")"
  printf 'artifact_wasm_sha256=%s\n' "$(sha256_file "$PUBLIC_DIR/$OUTPUT_BASENAME.wasm")"
  while IFS= read -r line; do
    printf 'fmm_%s\n' "$line"
  done < "$FMM3D_STAMP"
} > "$provenance_next"
mv -f "$provenance_next" "$BUILD_DIR/provenance.stamp"
echo "WASM_BUILD_OK architecture=$architecture lines=$lines bytes=$(wc -c < "$PUBLIC_DIR/$OUTPUT_BASENAME.wasm" | tr -d ' ')"
