/*
 * Emscripten bridge for the Fortran solver progress callback.
 *
 * The real scalar Fortran core (compiled under BIESOLVER_WASM_PROGRESS) emits
 * compact numeric progress events through a single bind(C) symbol,
 * `biesolver_progress`.  This translation unit resolves that symbol at the
 * final Emscripten link with an EM_JS definition that forwards each event to
 * `Module.onSolverProgress(...)`.
 *
 * Progress reporting is strictly observational: a missing or throwing
 * JavaScript handler is swallowed so it can never abort or perturb the
 * numerical solve.  This file carries no formatting, no solver logic, and no
 * mutable global state.
 */
#include <emscripten.h>
#include <stdint.h>

/*
 * With -sWASM_BIGINT=1 the four i64 payloads arrive in JavaScript as BigInt.
 * The handler is looked up on Module for every event so the Worker can install
 * and clear it per solve.  Any lookup or handler failure is swallowed.
 */
EM_JS(void, biesolver_progress,
      (int32_t stage, int64_t current, int64_t total,
       int64_t aux0, int64_t aux1, double value), {
  try {
    const handler = Module['onSolverProgress'];
    if (typeof handler === 'function') {
      handler(stage, current, total, aux0, aux1, value);
    }
  } catch (_) {
    // Progress reporting is observational and must never abort the solver.
  }
});
