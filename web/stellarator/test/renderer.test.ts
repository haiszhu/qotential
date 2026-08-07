import { describe, expect, it } from 'vitest';
import { prepareGpuMesh } from '../src/renderer';

describe('prepareGpuMesh', () => {
  it('centers and uniformly scales positions while preserving topology', () => {
    const mesh = prepareGpuMesh(
      new Float32Array([1, 2, 3, 5, 4, 3]),
      new Float32Array([-4, -2]),
      new Uint32Array([0, 1, 0]),
    );
    expect(mesh.positions).toHaveLength(6);
    expect(Math.max(...mesh.positions.map(Math.abs))).toBeCloseTo(1, 6);
    expect([...mesh.triangles]).toEqual([0, 1, 0]);
    expect(mesh.colors).toHaveLength(8);
  });
});
