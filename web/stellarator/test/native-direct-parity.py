#!/usr/bin/env python3
"""Validate deterministic native direct-solver fixtures."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import shutil
import struct
import subprocess
import tempfile
from pathlib import Path


WEB_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_BINARY = WEB_ROOT / "build-native" / "stellarator_direct"
FIXTURE_DIR = WEB_ROOT / "fixtures" / "native"
MAGIC = b"STGRF001"
FIELDS = [
    "sx", "snx", "sw", "ub", "ubn", "u",
    "render_xyz", "render_log_error", "render_triangles",
]


def run_fixture(binary: Path, output: Path, mp: int, np: int, order: int,
                surface: str, restol: float) -> str:
    completed = subprocess.run(
        [str(binary), str(output), str(mp), str(np), str(order),
         surface, repr(restol)],
        check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    if not output.is_file() or output.stat().st_size == 0:
        raise SystemExit(f"{binary}: did not create {output}")
    return completed.stdout


def parse_fixture(
    path: Path, mp: int, np: int, order: int, surface: str, restol: float,
) -> tuple[dict[str, object], dict[str, list[float] | list[int]]]:
    blob = path.read_bytes()
    if len(blob) < 48:
        raise SystemExit(f"{path}: truncated header")
    magic, ntri, nsrc, nrender, nfaces, grf_error = struct.unpack_from("<8s4qd", blob)
    if magic != MAGIC:
        raise SystemExit(f"{path}: bad magic {magic!r}")
    actual = {"ntri": ntri, "nsrc": nsrc, "nrender": nrender, "nfaces": nfaces}

    shapes = {
        "sx": 3 * nsrc, "snx": 3 * nsrc, "sw": nsrc,
        "ub": nsrc, "ubn": nsrc, "u": nsrc,
        "render_xyz": 3 * nrender, "render_log_error": nrender,
        "render_triangles": 3 * nfaces,
    }
    offset = 48
    arrays: dict[str, list[float] | list[int]] = {}
    offsets: dict[str, dict[str, int]] = {}
    for name in FIELDS:
        count = shapes[name]
        size = 8 * count
        if offset + size > len(blob):
            raise SystemExit(f"{path}: truncated {name}")
        code = "q" if name == "render_triangles" else "d"
        arrays[name] = list(struct.unpack_from(f"<{count}{code}", blob, offset))
        offsets[name] = {"offset": offset, "bytes": size, "count": count}
        offset += size
    if offset != len(blob):
        raise SystemExit(f"{path}: {len(blob) - offset} trailing bytes")

    for name, values in arrays.items():
        if name != "render_triangles" and not all(math.isfinite(value) for value in values):
            raise SystemExit(f"{path}: non-finite value in {name}")
    triangles = arrays["render_triangles"]
    if not all(0 <= value < nrender for value in triangles):
        raise SystemExit(f"{path}: triangle index out of range")
    logs = arrays["render_log_error"]
    if min(logs) < -16.0:
        raise SystemExit(f"{path}: render log floor below -16")
    measured = max(abs(value) for value in arrays["u"]) / max(
        abs(value) for value in arrays["ubn"]
    )
    if measured != grf_error:
        raise SystemExit(f"{path}: GRF scalar {grf_error} != measured {measured}")

    metadata: dict[str, object] = {
        "schema": 1, "mp": mp, "np": np, "order": order,
        "surface": surface, "restol": restol, "isimd": 0,
        "targetBlock": 64, "float": "little-endian-f64",
        "integer": "little-endian-i64", "fields": FIELDS,
        **actual, "grfError": grf_error, "offsets": offsets,
        "bytes": len(blob), "sha256": hashlib.sha256(blob).hexdigest(),
        "tolerances": {"nativeRepeatAbsolute": 0.0},
    }
    return metadata, arrays


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", type=Path, default=DEFAULT_BINARY)
    parser.add_argument("--mp", type=int, default=12)
    parser.add_argument("--np", type=int, default=36)
    parser.add_argument("--order", type=int, choices=range(4, 17, 2), default=4)
    parser.add_argument("--surface", choices=("builtin", "w7x"), default="builtin")
    parser.add_argument("--restol", type=float, default=1.0e-1)
    parser.add_argument("--update", action="store_true")
    args = parser.parse_args()
    if args.mp <= 0 or args.np <= 0:
        parser.error("--mp and --np must be positive integers")
    if args.restol <= 0.0:
        parser.error("--restol must be positive")
    if not args.binary.is_file():
        raise SystemExit(f"missing native direct solver: {args.binary}")
    stem = (f"order{args.order}-direct" if args.surface == "builtin"
            else f"w7x-order{args.order}-direct")
    fixture_bin = FIXTURE_DIR / f"{stem}.bin"
    fixture_json = FIXTURE_DIR / f"{stem}.json"

    with tempfile.TemporaryDirectory(prefix="stellarator-direct-") as tmp:
        first, second = Path(tmp) / "first.bin", Path(tmp) / "second.bin"
        first_log = run_fixture(args.binary, first, args.mp, args.np,
                                args.order, args.surface, args.restol)
        run_fixture(args.binary, second, args.mp, args.np, args.order,
                    args.surface, args.restol)
        if first.read_bytes() != second.read_bytes():
            raise SystemExit("two native direct runs are not bit-for-bit identical")
        metadata, _ = parse_fixture(first, args.mp, args.np, args.order,
                                    args.surface, args.restol)
        metadata["nativeOutput"] = first_log.strip()
        compiler = os.environ.get("FC", "gfortran-16")
        metadata["compilerVersion"] = subprocess.run(
            [compiler, "--version"], check=True, text=True,
            stdout=subprocess.PIPE,
        ).stdout.splitlines()[0]
        if args.update:
            FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(first, fixture_bin)
            fixture_json.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        else:
            if not fixture_bin.is_file() or first.read_bytes() != fixture_bin.read_bytes():
                raise SystemExit(f"native result differs from frozen {fixture_bin.name}")

            frozen = json.loads(fixture_json.read_text())
            if frozen["sha256"] != metadata["sha256"]:
                raise SystemExit("native metadata SHA differs from frozen fixture")
            for name in ("mp", "np", "order", "surface", "restol",
                         "ntri", "nsrc", "nrender", "nfaces"):
                if frozen.get(name) != metadata[name]:
                    raise SystemExit(
                        f"native metadata {name}={metadata[name]} differs from "
                        f"frozen value {frozen.get(name)}"
                    )
    print(
        "NATIVE_DIRECT_OK",
        f"ntri={metadata['ntri']}", f"nsrc={metadata['nsrc']}",
        f"grf={metadata['grfError']:.17e}", f"sha256={metadata['sha256']}",
    )


if __name__ == "__main__":
    main()
