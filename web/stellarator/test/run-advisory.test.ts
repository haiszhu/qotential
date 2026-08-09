import { describe, expect, it } from 'vitest';

import { MEASURED_W7X_REFERENCE, runAdvisory } from '../src/run-advisory';

describe('runAdvisory', () => {
  it('exports the measured W7-X reference', () => {
    expect(MEASURED_W7X_REFERENCE).toEqual({
      mp: 12,
      np: 36,
      order: 6,
      restol: 0.1,
      nsrc: 115_500,
    });
  });

  it('has no advisory for the built-in surface', () => {
    expect(runAdvisory('builtin', 4, 12, 36, 0.1, 'fmm', 1e-3)).toBe('');
  });

  it('warns that Direct is quadratic and ignores the FMM tolerance', () => {
    expect(runAdvisory('builtin', 4, 12, 36, 0.1, 'direct', 1e-9)).toMatch(
      /Direct.*O\(N²\).*1e-9.*ignored/i,
    );
  });

  it('identifies the missing order-4 W7-X parity reference', () => {
    expect(runAdvisory('w7x', 4, 12, 36, 0.1, 'fmm', 1e-3)).toMatch(
      /order 4.*no frozen W7-X parity reference.*order.?6.*measured reference/i,
    );
  });

  it('reports the exact measured W7-X source count', () => {
    expect(runAdvisory('w7x', 6, 12, 36, 0.1, 'fmm', 1e-3)).toMatch(
      /measured W7-X reference.*N=115,500.*115,500 source nodes/i,
    );
  });

  it('does not extrapolate unmeasured W7-X source counts', () => {
    expect(runAdvisory('w7x', 6, 12, 36, 0.03, 'fmm', 1e-3)).toMatch(
      /node count is data-dependent.*known after real geometry/i,
    );
  });

  it('composes the Direct warning with the W7-X reference advisory', () => {
    const text = runAdvisory('w7x', 6, 12, 36, 0.1, 'direct', 1e-6);
    expect(text).toMatch(/measured W7-X reference.*N=115,500/i);
    expect(text).toMatch(/Direct.*O\(N²\).*1e-6.*ignored/i);
  });
});
