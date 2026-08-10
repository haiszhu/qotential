import { describe, expect, it, vi } from 'vitest';
import {
  createSolverRuntimeSelector,
  detectMemory64,
  selectSolverRuntime,
} from '../src/memory64';

describe('Memory64 runtime selection', () => {
  it('selects wasm64 only after validation and instantiation succeed', async () => {
    const validate = vi.fn(() => true);
    const instantiate = vi.fn(async () => ({}));
    await expect(detectMemory64({ validate, instantiate })).resolves.toBe(true);
    expect(validate).toHaveBeenCalledOnce();
    expect(instantiate).toHaveBeenCalledOnce();
    expect(await selectSolverRuntime(
      { wasm32ModuleUrl: 'solver.js', wasm64ModuleUrl: 'solver64.js' },
      async () => true,
    )).toEqual({
      moduleUrl: 'solver64.js',
      pointerBits: 64,
      logLine: '[runtime] WebAssembly Memory64 enabled',
    });
  });

  it('falls back to wasm32 only when validation reports unsupported', async () => {
    const instantiate = vi.fn(async () => ({}));
    await expect(detectMemory64({ validate: () => false, instantiate }))
      .resolves.toBe(false);
    expect(instantiate).not.toHaveBeenCalled();
    expect(await selectSolverRuntime(
      { wasm32ModuleUrl: 'solver.js', wasm64ModuleUrl: 'solver64.js' },
      async () => false,
    )).toEqual({
      moduleUrl: 'solver.js',
      pointerBits: 32,
      logLine: '[runtime] Memory64 unavailable; using wasm32',
    });
  });

  it('does not disguise an unexpected instantiation failure as unsupported', async () => {
    const failure = new Error('probe instantiation failed');
    await expect(detectMemory64({
      validate: () => true,
      instantiate: async () => { throw failure; },
    })).rejects.toBe(failure);
  });

  it('probes once per selector and a fresh selector probes again', async () => {
    const probe = vi.fn(async () => true);
    const urls = { wasm32ModuleUrl: 'solver.js', wasm64ModuleUrl: 'solver64.js' };
    const firstWorker = createSolverRuntimeSelector(probe);
    await firstWorker(urls);
    await firstWorker(urls);
    expect(probe).toHaveBeenCalledOnce();

    const replacementWorker = createSolverRuntimeSelector(probe);
    await replacementWorker(urls);
    expect(probe).toHaveBeenCalledTimes(2);
  });
});
