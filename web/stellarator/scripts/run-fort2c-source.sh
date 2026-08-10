#!/usr/bin/env bash
set -euo pipefail

FORT2C_TOOL=${FORT2C_TOOL:-fort2c}
FORT2C_SRC=${FORT2C_SRC:-$HOME/fort2c}
command -v "$FORT2C_TOOL" >/dev/null 2>&1 || {
  echo "required tool not found: $FORT2C_TOOL" >&2
  exit 1
}
test -f "$FORT2C_SRC/fort2c/transpiler.py" || {
  echo "missing fort2c source package under $FORT2C_SRC" >&2
  exit 1
}

tool_path=$(command -v "$FORT2C_TOOL")
tool_real=$(python3 - "$tool_path" <<'PY'
import os
import sys
print(os.path.realpath(sys.argv[1]))
PY
)
tool_python=$(sed -n '1s/^#!//p' "$tool_real")
test -x "$tool_python" || {
  echo "cannot resolve fort2c Python interpreter from $tool_real" >&2
  exit 1
}

export PYTHONPATH="$FORT2C_SRC${PYTHONPATH:+:$PYTHONPATH}"
exec "$tool_python" -m fort2c.cli "$@"
