#!/usr/bin/env python3
"""Generate the Fortran W7-X mode module from the production MATLAB table."""

from __future__ import annotations

import argparse
import hashlib
import re
import struct
from pathlib import Path


HERE = Path(__file__).resolve()
WEB_ROOT = HERE.parents[1]
QOTENTIAL_ROOT = HERE.parents[3]
SOURCE = QOTENTIAL_ROOT / "test" / "w7x_grf.m"
MODULE = WEB_ROOT / "fortran" / "w7x_modes_mod.f90"
FIXTURE = WEB_ROOT / "geometry" / "fixtures" / "w7x.bin"


def parse_modes(path: Path) -> list[tuple[int, int, str, str]]:
    source = path.read_text()
    match = re.search(r"\nD\s*=\s*\[(.*?)\];", source, re.S)
    if match is None:
        raise SystemExit(f"{path}: no `D = [ ... ];` mode table found")
    body = re.sub(r"%[^\n]*", "", re.sub(r"\.\.\.\s*\n", "", match.group(1)))
    rows = [row.strip() for row in body.replace("\n", ";").split(";") if row.strip()]
    modes: list[tuple[int, int, str, str]] = []
    for row in rows:
        fields = [field for field in re.split(r"[,\s]+", row) if field]
        if len(fields) != 4:
            raise SystemExit(f"{path}: expected four columns, got {len(fields)} in {row!r}")
        modes.append((int(fields[0]), int(fields[1]), fields[2], fields[3]))
    if len(modes) != 288:
        raise SystemExit(f"{path}: expected 288 modes, got {len(modes)}")
    if float(modes[0][2]) != 5.51984749975020073:
        raise SystemExit(f"{path}: unexpected R00 {modes[0][2]}")
    return modes


def fixture_bytes(modes: list[tuple[int, int, str, str]]) -> bytes:
    nmode = len(modes)
    columns = (
        [float(mode[0]) for mode in modes],
        [float(mode[1]) for mode in modes],
        [float(mode[2]) for mode in modes],
        [float(mode[3]) for mode in modes],
    )
    return struct.pack("<q", nmode) + b"".join(
        struct.pack(f"<{nmode}d", *column) for column in columns
    )


def wrapped(values: list[str], per_line: int) -> str:
    lines = []
    for offset in range(0, len(values), per_line):
        chunk = ", ".join(values[offset : offset + per_line])
        suffix = ", &" if offset + per_line < len(values) else " ]"
        lines.append(f"    {chunk}{suffix}")
    return "\n".join(lines)


def scalar_assignments(name: str, values: list[str], kind: str) -> str:
    return "\n".join(
        f"    {name}({index}) = {value}_{kind}"
        for index, value in enumerate(values, start=1)
    )


def module_text(modes: list[tuple[int, int, str, str]]) -> str:
    mn = [f"{value}_8" for mode in modes for value in mode[:2]]
    rc = [f"{mode[2]}_r64" for mode in modes]
    zs = [f"{mode[3]}_r64" for mode in modes]
    mn_plain = [str(value) for mode in modes for value in mode[:2]]
    rc_plain = [mode[2] for mode in modes]
    zs_plain = [mode[3] for mode in modes]
    return f"""! Generated from test/w7x_grf.m by generate-w7x-modes.py.
! Do not edit this coefficient table by hand.
module w7x_modes_mod
  use quatapproximation_mod, only: r64
  implicit none
  private
  integer(8), parameter, public :: W7X_NMODE = 288_8
  integer(8), parameter, public :: W7X_NFP = 5_8
  public :: load_w7x_modes
  integer(8), parameter, public :: W7X_MN(2*W7X_NMODE) = [ &
{wrapped(mn, 8)}
  real(r64), parameter, public :: W7X_RC(W7X_NMODE) = [ &
{wrapped(rc, 3)}
  real(r64), parameter, public :: W7X_ZS(W7X_NMODE) = [ &
{wrapped(zs, 3)}
contains
  subroutine load_w7x_modes(mn, rc, zs)
    real(r64), intent(out) :: mn(2*W7X_NMODE), rc(W7X_NMODE), zs(W7X_NMODE)
{scalar_assignments('mn', mn_plain, 'r64')}
{scalar_assignments('rc', rc_plain, 'r64')}
{scalar_assignments('zs', zs_plain, 'r64')}
  end subroutine load_w7x_modes
end module w7x_modes_mod
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="verify generated files are current")
    args = parser.parse_args()

    modes = parse_modes(SOURCE)
    expected_module = module_text(modes)
    expected_fixture = fixture_bytes(modes)
    if args.check:
        if not MODULE.is_file() or MODULE.read_text() != expected_module:
            raise SystemExit(f"{MODULE}: missing or stale; run {HERE}")
        if not FIXTURE.is_file() or FIXTURE.read_bytes() != expected_fixture:
            raise SystemExit(f"{FIXTURE}: does not match {SOURCE}")
    else:
        MODULE.write_text(expected_module)

    digest = hashlib.sha256(expected_fixture).hexdigest()
    print(f"W7X_MODES_OK nmode={len(modes)} sha256={digest}")


if __name__ == "__main__":
    main()
