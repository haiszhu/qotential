#!/usr/bin/env python3
"""Restore subroutine-local array lifetime in LFortran C output.

LFortran's C backend currently emits plain malloc calls for automatic arrays
inside Fortran subroutines, but does not emit matching frees when the
subroutine returns. In a long-lived WASM module those allocations accumulate
until wasm32's address space is exhausted.

For void C functions (the representation used for Fortran subroutines), move
those local allocations to the solver's resettable scratch arena. Deliberately
leave non-void functions alone: generated array-return helpers may return
heap-backed data whose lifetime extends beyond the helper call.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path


VOID_DEFINITION = re.compile(r"^void\s+[A-Za-z_]\w*\s*\(")
COMPILER_TEMP = re.compile(r"\b(__(?:libasr|lcompilers)[A-Za-z0-9_]*)->data\b")
RUNTIME_ALLOC = re.compile(
    r"_lfortran_malloc_alloc\(_lfortran_get_default_allocator\(\), (.*)\)"
)
RUNTIME_FREE = re.compile(
    r"_lfortran_free_alloc\(_lfortran_get_default_allocator\(\), "
    r"\(char\*\)(__(?:libasr|lcompilers)[A-Za-z0-9_]*)->data\)"
)

PATCH_LEVEL_CLEAR = "void clear_patch_levels_r64(struct patch_levels_t* lv)"
PATCH_LEVEL_STORAGE_FREE = "    if (lv->lev->is_allocated) {"
PATCH_LEVEL_DESCRIPTOR_CLEANUP = """\
    /* LFortran represents each nested allocatable component by a separately
       allocated descriptor. Fortran DEALLOCATE releases only descriptor->data
       in the current C backend, so release the five descriptors explicitly. */
    for (int64_t __biesolver_level = 0;
         __biesolver_level < lv->lev->dims[0].length; ++__biesolver_level) {
        _lfortran_free_alloc(_lfortran_get_default_allocator(),
            (char*)lv->lev->data[__biesolver_level].sx);
        _lfortran_free_alloc(_lfortran_get_default_allocator(),
            (char*)lv->lev->data[__biesolver_level].snx);
        _lfortran_free_alloc(_lfortran_get_default_allocator(),
            (char*)lv->lev->data[__biesolver_level].sw);
        _lfortran_free_alloc(_lfortran_get_default_allocator(),
            (char*)lv->lev->data[__biesolver_level].qpoint);
        _lfortran_free_alloc(_lfortran_get_default_allocator(),
            (char*)lv->lev->data[__biesolver_level].qradii2);
    }
"""


def rewrite(source: str) -> tuple[str, int]:
    lines = source.splitlines(keepends=True)
    output: list[str] = []
    pending_void = False
    in_void = False
    brace_depth = 0
    replacements = 0

    for line in lines:
        before: list[str] = []
        after: list[str] = []
        stripped = line.strip()
        if not in_void and VOID_DEFINITION.match(line) and not stripped.endswith(";"):
            pending_void = True

        if pending_void or in_void:
            count = line.count("malloc(")
            if count:
                line = line.replace("malloc(", "biesolver_scope_alloc(")
                replacements += count

            opens = line.count("{")
            closes = line.count("}")
            if pending_void and opens:
                in_void = True
                pending_void = False
                after.append(
                    "    struct biesolver_scope_mark __biesolver_scope = "
                    "biesolver_scope_enter();\n"
                )
            if in_void:
                temporary = COMPILER_TEMP.search(line)
                if temporary and "_lfortran_malloc_alloc(" in line:
                    line, count = RUNTIME_ALLOC.subn(
                        lambda match: f"biesolver_scope_alloc({match.group(1)})", line
                    )
                    if count:
                        replacements += count
                free = RUNTIME_FREE.search(line)
                if free:
                    line = RUNTIME_FREE.sub("(void)0", line)
                    replacements += 1
                if stripped == "return;":
                    before.append("    biesolver_scope_leave(__biesolver_scope);\n")
                next_depth = brace_depth + opens - closes
                if next_depth == 0:
                    final_close = line.rfind("}")
                    if final_close < 0:
                        raise ValueError("void function ended without a closing brace")
                    leading_closures = line[:final_close].rstrip()
                    if leading_closures:
                        before.append(leading_closures + "\n")
                    before.append("    biesolver_scope_leave(__biesolver_scope);\n")
                    line = line[final_close:]
                brace_depth = next_depth
                if brace_depth == 0:
                    in_void = False

        output.extend(before)
        output.append(line)
        output.extend(after)

    if pending_void or in_void or brace_depth != 0:
        raise ValueError("could not find the end of a generated void function")
    rewritten = "".join(output)
    # The generated file contains a prototype and then the definition.
    clear_start = rewritten.rfind(PATCH_LEVEL_CLEAR)
    if clear_start >= 0:
        storage_free = rewritten.find(PATCH_LEVEL_STORAGE_FREE, clear_start)
        if storage_free < 0:
            raise ValueError("could not find patch-level descriptor cleanup point")
        rewritten = (
            rewritten[:storage_free]
            + PATCH_LEVEL_DESCRIPTOR_CLEANUP
            + rewritten[storage_free:]
        )
        replacements += 1
    return rewritten, replacements


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=Path)
    args = parser.parse_args()

    original = args.path.read_text()
    rewritten, replacements = rewrite(original)
    if replacements == 0:
        raise SystemExit("no subroutine-local malloc calls found in LFortran C output")
    args.path.write_text(rewritten)
    print(f"LFORTRAN_C_LIFETIME_REWRITES={replacements}")


if __name__ == "__main__":
    main()
