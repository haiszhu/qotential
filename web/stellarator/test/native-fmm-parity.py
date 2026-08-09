#!/usr/bin/env python3
"""Validate temporary native FMM fixtures before WASM comparison."""

from __future__ import annotations

import argparse
import math
import struct
from pathlib import Path


MAGIC = b"STGRF001"
CASES = (
    "builtin-order4-e3.bin",
    "builtin-order6-e6.bin",
    "w7x-order6-e6.bin",
)


def validate(path: Path) -> tuple[int, int, int, int, float]:
    blob = path.read_bytes()
    if len(blob) < 48:
        raise SystemExit(f"{path}: truncated header")
    magic, ntri, nsrc, nrender, nfaces, grf = struct.unpack_from("<8s4qd", blob)
    if magic != MAGIC:
        raise SystemExit(f"{path}: bad magic {magic!r}")
    if min(ntri, nsrc, nrender, nfaces) <= 0 or not math.isfinite(grf):
        raise SystemExit(f"{path}: invalid header values")
    counts = (3 * nsrc, 3 * nsrc, nsrc, nsrc, nsrc, nsrc,
              3 * nrender, nrender)
    offset = 48
    for count in counts:
        end = offset + 8 * count
        if end > len(blob):
            raise SystemExit(f"{path}: truncated floating-point field")
        for (value,) in struct.iter_unpack("<d", memoryview(blob)[offset:end]):
            if not math.isfinite(value):
                raise SystemExit(f"{path}: non-finite floating-point field")
        offset = end
    triangle_count = 3 * nfaces
    end = offset + 8 * triangle_count
    if end != len(blob):
        raise SystemExit(f"{path}: fixture length mismatch")
    for (value,) in struct.iter_unpack("<q", memoryview(blob)[offset:end]):
        if value < 0 or value >= nrender:
            raise SystemExit(f"{path}: triangle index {value} out of range")
    return ntri, nsrc, nrender, nfaces, grf


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture-dir", type=Path, required=True)
    args = parser.parse_args()
    for name in CASES:
        path = args.fixture_dir / name
        if not path.is_file():
            raise SystemExit(f"missing native FMM fixture: {path}")
        ntri, nsrc, nrender, nfaces, grf = validate(path)
        print("NATIVE_FMM_FIXTURE_OK", name, f"ntri={ntri}", f"nsrc={nsrc}",
              f"nrender={nrender}", f"nfaces={nfaces}", f"grf={grf:.17e}")


if __name__ == "__main__":
    main()
