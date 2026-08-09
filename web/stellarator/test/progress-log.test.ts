import { describe, expect, it } from 'vitest';
import {
  createProgressForwarder,
  formatProgressEvent,
  safeCount,
  withSolverProgress,
  type RawProgressCallback,
} from '../src/progress-log';
import type { SolverWasmModule } from '../src/wasm-loader';

const NO_VALUE = 0;

describe('safeCount', () => {
  it('accepts counts within the safe-integer range', () => {
    expect(safeCount(0n, 'x')).toBe(0);
    expect(safeCount(115_500n, 'x')).toBe(115_500);
  });

  it('rejects counts outside Number.MAX_SAFE_INTEGER', () => {
    const tooBig = BigInt(Number.MAX_SAFE_INTEGER) + 1n;
    expect(() => safeCount(tooBig, 'nsrc')).toThrow(/safe-integer/i);
    expect(() => safeCount(-1n, 'nsrc')).toThrow(/safe-integer/i);
  });
});

describe('formatProgressEvent', () => {
  it('formats geometry begin differently for built-in and W7-X', () => {
    const builtin = formatProgressEvent(1, 0n, 0n, 0n, 0n, NO_VALUE)!;
    const w7x = formatProgressEvent(1, 0n, 0n, 1n, 0n, NO_VALUE)!;
    expect(builtin.line).not.toBe(w7x.line);
    expect(w7x.line).toMatch(/W7-X adaptive charts/);
    expect(builtin.alwaysForward).toBe(true);
    expect(w7x.alwaysForward).toBe(true);
  });

  it('formats geometry ready with panel and source-node counts', () => {
    const event = formatProgressEvent(2, 5_500n, 115_500n, 4n, 6n, NO_VALUE)!;
    expect(event.line).toBe('[geometry] 5,500 panels, 115,500 source nodes');
    expect(event.stageKey).toBe('geometry');
  });

  it('formats a direct milestone with right-aligned percent', () => {
    const event = formatProgressEvent(4, 37n, 1_805n, 0n, 0n, NO_VALUE)!;
    expect(event.line).toBe('[direct]   2%  block 37 / 1,805');
    expect(event.alwaysForward).toBe(false);
  });

  it('formats FMM begin and one real-time completion event', () => {
    const begin = formatProgressEvent(3, 0n, 1n, 1n, 0n, 0)!;
    const complete = formatProgressEvent(4, 1n, 1n, 1n, 0n, 0.125)!;
    expect(begin.line).toMatch(/\[fmm\].*Laplace GRF/i);
    expect(complete.line).toMatch(/\[fmm\].*completed.*0\.125 s/i);
    expect(begin.stageKey).toBe('fmm');
    expect(complete.stageKey).toBe('fmm');
    expect(begin.alwaysForward).toBe(true);
    expect(complete.alwaysForward).toBe(true);
  });

  it('ignores unknown far-field kernel payloads', () => {
    expect(formatProgressEvent(3, 0n, 1n, 2n, 0n, 0)).toBeUndefined();
    expect(formatProgressEvent(4, 1n, 1n, -1n, 0n, 0)).toBeUndefined();
  });

  it('marks a 100-percent close-count milestone as always forwarded', () => {
    const event = formatProgressEvent(6, 5_500n, 5_500n, 0n, 0n, NO_VALUE)!;
    expect(event.line).toBe('[close/count] 100%  patch 5,500 / 5,500');
    expect(event.alwaysForward).toBe(true);
  });

  it('shares one stageKey between a phase begin and its progress', () => {
    expect(formatProgressEvent(3, 0n, 1_805n, 0n, 0n, NO_VALUE)!.stageKey).toBe('direct');
    expect(formatProgressEvent(4, 37n, 1_805n, 0n, 0n, NO_VALUE)!.stageKey).toBe('direct');
    expect(formatProgressEvent(7, 0n, 5_500n, 0n, 0n, NO_VALUE)!.stageKey).toBe('close/rrq');
    expect(formatProgressEvent(8, 110n, 5_500n, 0n, 0n, NO_VALUE)!.stageKey).toBe('close/rrq');
  });

  it('formats the result line with a signed two-digit exponent', () => {
    const event = formatProgressEvent(11, 44_280n, 20_301n, 0n, 0n, 2.31594e-3)!;
    expect(event.line).toBe('[result] GRF max relative error = 2.315940e-03');
    expect(event.alwaysForward).toBe(true);
  });

  it('ignores unknown stage codes', () => {
    expect(formatProgressEvent(0, 0n, 0n, 0n, 0n, NO_VALUE)).toBeUndefined();
    expect(formatProgressEvent(99, 0n, 0n, 0n, 0n, NO_VALUE)).toBeUndefined();
  });

  it('propagates a safe-integer violation from the decoder', () => {
    const tooBig = BigInt(Number.MAX_SAFE_INTEGER) + 1n;
    expect(() => formatProgressEvent(2, tooBig, 1n, 0n, 0n, NO_VALUE)).toThrow(/safe-integer/i);
  });
});

describe('createProgressForwarder', () => {
  function harness() {
    let clock = 0;
    const lines: string[] = [];
    const forward = createProgressForwarder({ sink: (line) => lines.push(line), now: () => clock });
    return {
      lines,
      at(ms: number) { clock = ms; },
      forward,
    };
  }

  it('throttles ordinary milestones to one per 250 ms per stage', () => {
    const h = harness();
    // Begin line records the phase timestamp at t=0.
    h.at(0);
    h.forward(3, 0n, 1_805n, 0n, 0n, NO_VALUE);
    // Ordinary milestones inside the 250 ms window are suppressed.
    h.at(0);
    h.forward(4, 37n, 1_805n, 0n, 0n, NO_VALUE);
    h.at(100);
    h.forward(4, 73n, 1_805n, 0n, 0n, NO_VALUE);
    // The next ordinary milestone at 250 ms is forwarded.
    h.at(250);
    h.forward(4, 109n, 1_805n, 0n, 0n, NO_VALUE);
    expect(h.lines).toEqual([
      '[direct] evaluating naive Laplace GRF',
      '[direct]   6%  block 109 / 1,805',
    ]);
  });

  it('never throttles begin, 100-percent, scatter, render, and result', () => {
    const h = harness();
    h.at(0);
    h.forward(5, 0n, 5_500n, 0n, 0n, NO_VALUE);       // begin
    h.at(1);
    h.forward(6, 5_500n, 5_500n, 0n, 0n, NO_VALUE);   // 100 percent
    h.at(2);
    h.forward(9, 5_500n, 5_500n, 0n, 0n, NO_VALUE);   // scatter
    h.at(3);
    h.forward(10, 5_500n, 0n, 0n, 0n, NO_VALUE);      // render
    h.at(4);
    h.forward(11, 44_280n, 20_301n, 0n, 0n, 1.5e-2);  // result
    expect(h.lines).toEqual([
      '[close/count] finding close targets',
      '[close/count] 100%  patch 5,500 / 5,500',
      '[close/scatter] accumulating corrected values',
      '[render] building visualization lattice',
      '[result] GRF max relative error = 1.500000e-02',
    ]);
  });

  it('keeps throttle state independent per stage', () => {
    const h = harness();
    h.at(0);
    h.forward(4, 37n, 1_805n, 0n, 0n, NO_VALUE);      // direct ordinary, forwarded
    h.at(10);
    h.forward(6, 110n, 5_500n, 0n, 0n, NO_VALUE);     // close/count ordinary, forwarded (own stage)
    h.at(20);
    h.forward(4, 73n, 1_805n, 0n, 0n, NO_VALUE);      // direct within 250 ms, suppressed
    expect(h.lines).toEqual([
      '[direct]   2%  block 37 / 1,805',
      '[close/count]   2%  patch 110 / 5,500',
    ]);
  });

  it('ignores unknown stage codes and never throws on bad payloads', () => {
    const h = harness();
    h.forward(0, 0n, 0n, 0n, 0n, NO_VALUE);           // unknown -> ignored
    const tooBig = BigInt(Number.MAX_SAFE_INTEGER) + 1n;
    expect(() => h.forward(2, tooBig, 1n, 0n, 0n, NO_VALUE)).not.toThrow(); // decode error swallowed
    expect(h.lines).toEqual([]);
  });
});

describe('withSolverProgress', () => {
  function fakeModule(): SolverWasmModule {
    return {} as unknown as SolverWasmModule;
  }
  const noop: RawProgressCallback = () => {};

  it('clears the module callback after a successful run', () => {
    const wasm = fakeModule();
    const result = withSolverProgress(wasm, noop, () => {
      expect(wasm.onSolverProgress).toBe(noop);
      return 42;
    });
    expect(result).toBe(42);
    expect(wasm.onSolverProgress).toBeUndefined();
  });

  it('clears the module callback after a throwing run', () => {
    const wasm = fakeModule();
    expect(() =>
      withSolverProgress(wasm, noop, () => {
        throw new Error('solve failed');
      }),
    ).toThrow(/solve failed/);
    expect(wasm.onSolverProgress).toBeUndefined();
  });
});
