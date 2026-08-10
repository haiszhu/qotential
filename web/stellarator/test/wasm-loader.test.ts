import { describe, expect, it } from 'vitest';
import {
  checkedHeapOffset,
  createSimplexPrecomp,
  pointerArgument,
  requireHeapRange,
  requireStatus,
  SolverStatusError,
  type SolverWasmModule,
} from '../src/wasm-loader';

describe('WASM pointer normalization', () => {
  it('keeps wasm32 offsets numeric and converts wasm64 arguments to bigint', () => {
    expect(pointerArgument(64, 32)).toBe(64);
    expect(pointerArgument(64, 64)).toBe(64n);
    expect(() => pointerArgument(-1, 64)).toThrow(RangeError);
  });

  it('checks pointer results and heap view ranges before conversion', () => {
    expect(checkedHeapOffset(64)).toBe(64);
    expect(checkedHeapOffset(64n)).toBe(64);
    expect(() => checkedHeapOffset(BigInt(Number.MAX_SAFE_INTEGER) + 1n))
      .toThrow(RangeError);
    const heap = { HEAPU8: new Uint8Array(128) };
    expect(() => requireHeapRange(heap, 64, 64)).not.toThrow();
    expect(() => requireHeapRange(heap, 65, 64)).toThrow(RangeError);
  });
});

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

describe('simplex precomputation ownership', () => {
  function mockModule(options: { failAllocationAt?: number; status?: number } = {}) {
    const allocations: number[] = [];
    const frees: number[] = [];
    const calls: unknown[][] = [];
    let nextPointer = 64;
    const wasm = {
      HEAPU8: new Uint8Array(1 << 20),
      _malloc(bytes: number) {
        allocations.push(bytes);
        if (allocations.length === options.failAllocationAt) return 0;
        const pointer = nextPointer;
        nextPointer += bytes;
        return pointer;
      },
      _free(pointer: number) { frees.push(pointer); },
      _solver_simplex_precomp(...args: unknown[]) {
        calls.push(args);
        return options.status ?? 0;
      },
      _solver_last_error() { return options.status ?? 0; },
    } as unknown as SolverWasmModule;
    return { wasm, allocations, frees, calls };
  }

  it('allocates, initializes, and idempotently frees all seven arrays', () => {
    const state = mockModule();
    const refs = createSimplexPrecomp(state.wasm, 6);
    expect(state.allocations).toEqual([
      8 * 8, 8 * 8, 8 * 8 * 8, 8 * 8, 8 * 8 * 8,
      8 * 21 * 21, 8 * 21 * 21,
    ]);
    expect(state.calls).toHaveLength(1);
    expect(state.calls[0]?.slice(0, 3)).toEqual([8n, 5n, 21n]);
    refs.dispose();
    refs.dispose();
    expect(state.frees).toHaveLength(7);
  });

  it('passes every owned pointer through the wasm64 bigint ABI', () => {
    const state = mockModule();
    const refs = createSimplexPrecomp(state.wasm, 4, 64);
    expect(state.calls[0]?.slice(3).every((value) => typeof value === 'bigint')).toBe(true);
    refs.dispose();
    expect(state.frees.every((value) => typeof value === 'bigint')).toBe(true);
  });

  it('frees partial allocations and all arrays after precomputation failure', () => {
    const partial = mockModule({ failAllocationAt: 4 });
    expect(() => createSimplexPrecomp(partial.wasm, 6)).toThrow(/allocation/i);
    expect(partial.frees).toHaveLength(3);

    const failed = mockModule({ status: 109 });
    expect(() => createSimplexPrecomp(failed.wasm, 6)).toThrow(SolverStatusError);
    expect(failed.frees).toHaveLength(7);
  });
});
