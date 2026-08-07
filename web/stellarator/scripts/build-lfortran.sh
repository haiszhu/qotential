#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SOURCE_TREE=${1:-${LF_SRC:-}}
PINNED_COMMIT=bc9740b2d4a946e26418473262b19b06eb268427
PATCH="$SCRIPT_DIR/../patches/lfortran-c-backend.patch"

if [[ -z $SOURCE_TREE ]]; then
  echo "usage: $0 /path/to/lfortran-source" >&2
  echo "or set LF_SRC to that checkout" >&2
  exit 2
fi
test -d "$SOURCE_TREE/.git"
test -s "$PATCH"

actual=$(git -C "$SOURCE_TREE" rev-parse HEAD)
if [[ $actual != "$PINNED_COMMIT" ]]; then
  echo "LFortran checkout is $actual; expected $PINNED_COMMIT" >&2
  exit 1
fi

if git -C "$SOURCE_TREE" apply --reverse --check "$PATCH" >/dev/null 2>&1; then
  echo "LFortran C-backend patch is already applied"
elif git -C "$SOURCE_TREE" apply --check "$PATCH"; then
  git -C "$SOURCE_TREE" apply "$PATCH"
  echo "applied LFortran C-backend patch"
else
  echo "LFortran checkout is neither clean nor patched as expected" >&2
  exit 1
fi

unexpected=$(git -C "$SOURCE_TREE" diff --name-only | grep -Ev \
  '^(src/libasr/codegen/asr_to_c\.cpp|src/libasr/codegen/asr_to_c_cpp\.h|src/libasr/codegen/c_utils\.h|src/libasr/pass/intrinsic_subroutines\.h|src/libasr/runtime/lfortran_intrinsics\.h)$' || true)
if [[ -n $unexpected ]]; then
  echo "LFortran checkout has unrelated tracked changes:" >&2
  printf '%s\n' "$unexpected" >&2
  exit 1
fi

actual_patch=$(mktemp)
expected_patch=$(mktemp)
trap 'rm -f "$actual_patch" "$expected_patch"' EXIT
# Compare the applied change, not the rendering of it.  Two things differ
# between machines and neither is a difference in the code:
#   * blob abbreviation length in "index abc1234..def5678" -- git picks it from
#     core.abbrev (auto by default), which grows with the object count of the
#     particular clone;
#   * trailing whitespace on blank context lines -- a unified diff writes those
#     as a single space, which editors and pre-commit hooks strip out of the
#     committed patch file, after which git's own output can never match it.
normalise() { sed -E 's/^index [0-9a-f]+\.\.[0-9a-f]+/index HASH..HASH/; s/[[:space:]]+$//' "$1"; }
git -C "$SOURCE_TREE" diff --no-color -- \
  src/libasr/codegen/asr_to_c.cpp \
  src/libasr/codegen/asr_to_c_cpp.h \
  src/libasr/codegen/c_utils.h \
  src/libasr/pass/intrinsic_subroutines.h \
  src/libasr/runtime/lfortran_intrinsics.h | normalise /dev/stdin > "$actual_patch"
normalise "$PATCH" > "$expected_patch"
if ! cmp -s "$expected_patch" "$actual_patch"; then
  echo "applied compiler changes do not match the pinned patch" >&2
  diff -u "$expected_patch" "$actual_patch" | head -40 >&2
  exit 1
fi

if [[ ! -d $SOURCE_TREE/build ]]; then
  echo "no $SOURCE_TREE/build -- configure it first; see the macOS toolchain" >&2
  echo "section of web/stellarator/README.md for the cmake line" >&2
  exit 1
fi
# The default target, not --target lfortran: that one produces the binary
# without LFortran's intrinsic modfiles, and the probe below then fails with
# "Module 'lfortran_intrinsic_iso_c_binding' modfile was not found".
cmake --build "$SOURCE_TREE/build" -j "${JOBS:-4}"
test -x "$SOURCE_TREE/build/src/bin/lfortran"
"$SOURCE_TREE/build/src/bin/lfortran" --version
LF="$SOURCE_TREE/build/src/bin/lfortran" "$SCRIPT_DIR/../test/probe-script.test.sh"
