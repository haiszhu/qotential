#!/usr/bin/env bash
set -euo pipefail

WEB_ROOT=$(cd "$(dirname "$0")/.." && pwd)
FIXER="$WEB_ROOT/../../scripts/fix-lfortran-llvm.pl"

out=$(printf '%s\n' \
  'define void @probe() #0 {' \
  '  ret void' \
  '}' \
  'declare void @llvm.memset.p0.i64(ptr writeonly captures(none), i8, i64, i1 immarg)' \
  'attributes #0 = { nocallback nocreateundeforpoison nounwind }' | \
  perl "$FIXER")

grep -q 'ptr writeonly, i8' <<<"$out"
grep -q 'attributes #0 = { nocallback nounwind }' <<<"$out"
! grep -q 'nocreateundeforpoison' <<<"$out"
! grep -q 'captures(none)' <<<"$out"
