#!/usr/bin/env python3
"""Print a deterministic SHA-256 for Python translator source files."""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: source-tree-digest.py <package-directory>")
    root = Path(sys.argv[1]).resolve(strict=True)
    if not root.is_dir():
        raise SystemExit(f"not a directory: {root}")

    files = sorted(path for path in root.rglob("*.py") if "__pycache__" not in path.parts)
    if not files:
        raise SystemExit(f"no Python source under {root}")

    digest = hashlib.sha256()
    for path in files:
        resolved = path.resolve(strict=True)
        if os.path.commonpath((root, resolved)) != str(root):
            raise SystemExit(f"source escapes package tree: {path}")
        relative = path.relative_to(root).as_posix().encode("utf-8")
        payload = resolved.read_bytes()
        digest.update(len(relative).to_bytes(8, "little"))
        digest.update(relative)
        digest.update(len(payload).to_bytes(8, "little"))
        digest.update(payload)
    print(digest.hexdigest())


if __name__ == "__main__":
    main()
