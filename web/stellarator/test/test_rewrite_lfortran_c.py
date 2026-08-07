import importlib.util
import re
import unittest
from pathlib import Path


SCRIPT = Path(__file__).parents[1] / "scripts" / "rewrite-lfortran-c.py"
SPEC = importlib.util.spec_from_file_location("rewrite_lfortran_c", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class RewriteLFortranCTest(unittest.TestCase):
    def test_only_void_function_bodies_use_stack_lifetime(self):
        source = """\
void prototype(double *x);
void solver_step(double *x)
{
    double *tmp = (double*) malloc(sizeof(double) * 8);
    struct r64 temporary_value = {0};
    struct r64 *__libasr_created_temporary = &temporary_value;
    __libasr_created_temporary->data = (double*) _lfortran_malloc_alloc(_lfortran_get_default_allocator(), 8*sizeof(double));
    _lfortran_free_alloc(_lfortran_get_default_allocator(), (char*)__libasr_created_temporary->data);
    if (x == 0) {
        return;
    }
    if (x) {
        int *nested = (int*) malloc(sizeof(int) * 2);
    }
}
struct array *returning_helper(void)
{
    struct array *result = (struct array*) malloc(sizeof(struct array));
    return result;
}
"""
        rewritten, count = MODULE.rewrite(source)
        self.assertEqual(count, 4)
        self.assertIn("void prototype(double *x);", rewritten)
        self.assertEqual(rewritten.count("biesolver_scope_alloc("), 3)
        self.assertIn("(void)0;", rewritten)
        self.assertEqual(rewritten.count("biesolver_scope_enter()"), 1)
        self.assertEqual(rewritten.count("biesolver_scope_leave("), 2)
        self.assertIn("result = (struct array*) malloc(", rewritten)

    def test_nested_patch_descriptors_are_released(self):
        source = """\
void clear_patch_levels_r64(struct patch_levels_t* lv)
{
    if (lv->lev->is_allocated) {
        release_level_data(lv);
    }
}
"""
        rewritten, count = MODULE.rewrite(source)
        self.assertEqual(count, 1)
        self.assertEqual(rewritten.count("lv->lev->data[__biesolver_level]."), 5)
        self.assertLess(
            rewritten.index("__biesolver_level"),
            rewritten.index("if (lv->lev->is_allocated)"),
        )

    def test_scope_leave_is_outside_combined_closing_braces(self):
        source = """\
void combined_close(bool enabled)
{
    if (enabled) {
        use_scratch();
    }}
"""
        rewritten, _ = MODULE.rewrite(source)
        self.assertRegex(
            rewritten,
            re.compile(
                r"use_scratch\(\);\n\s*}\n"
                r"\s*biesolver_scope_leave\(__biesolver_scope\);\n}"
            ),
        )


if __name__ == "__main__":
    unittest.main()
