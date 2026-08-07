import { describe, expect, it } from 'vitest';
import { requireStatus, SolverStatusError } from '../src/wasm-loader';

describe('WASM status handling', () => {
  it('accepts zero and reports the solver error code for failures', () => {
    expect(() => requireStatus('prepare', 0, () => 0)).not.toThrow();
    expect(() => requireStatus('close', 7, () => 42)).toThrow(/close.*7.*42/i);
  });

  it('throws a structured SolverStatusError carrying phase, status, and solver error', () => {
    try {
      requireStatus('close', 7, () => 42);
      throw new Error('expected requireStatus to throw');
    } catch (error) {
      expect(error).toBeInstanceOf(SolverStatusError);
      const status = error as SolverStatusError;
      expect(status.phase).toBe('close');
      expect(status.status).toBe(7);
      expect(status.solverError).toBe(42);
      expect(status.message).toBe('close failed: status=7, solver error=42');
    }
  });
});
