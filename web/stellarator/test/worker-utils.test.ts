import { describe, expect, it } from 'vitest';
import { toSolverDataset } from '../src/worker-utils';

describe('toSolverDataset', () => {
  it('converts the WASM ABI mesh into browser buffers', () => {
    const data = toSolverDataset({
      nsrc: 2,
      nrender: 3,
      nfaces: 1,
      grfError: 1e-4,
      positions: new Float64Array([1, 2, 3, 4, 5, 6, 7, 8, 9]),
      logError: new Float64Array([-5, -4, -3]),
      triangles: new BigInt64Array([0n, 1n, 2n]),
    }, 17);

    expect(data.positions).toEqual(new Float32Array([1, 2, 3, 4, 5, 6, 7, 8, 9]));
    expect(data.logError).toEqual(new Float32Array([-5, -4, -3]));
    expect(data.triangles).toEqual(new Uint32Array([0, 1, 2]));
    expect(data.elapsedMs).toBe(17);
  });

  it('rejects an invalid Fortran mesh index', () => {
    expect(() => toSolverDataset({
      nsrc: 1,
      nrender: 2,
      nfaces: 1,
      grfError: 0,
      positions: new Float64Array(6),
      logError: new Float64Array(2),
      triangles: new BigInt64Array([0n, 1n, 2n]),
    }, 1)).toThrow(/triangle index/i);
  });
});
