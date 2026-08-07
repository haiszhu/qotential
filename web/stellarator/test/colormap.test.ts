import { describe, expect, it } from 'vitest';
import { colorizeLogError, LOG_ERROR_DOMAIN } from '../src/colormap';

describe('colorizeLogError', () => {
  it('uses the full solver floor-to-unity logarithmic domain', () => {
    expect(LOG_ERROR_DOMAIN).toEqual([-16, 0]);
  });

  it('clamps values to the fixed display domain', () => {
    const low = colorizeLogError(new Float32Array([-100]));
    const edge = colorizeLogError(new Float32Array([LOG_ERROR_DOMAIN[0]]));
    const high = colorizeLogError(new Float32Array([100]));
    const top = colorizeLogError(new Float32Array([LOG_ERROR_DOMAIN[1]]));
    expect([...low]).toEqual([...edge]);
    expect([...high]).toEqual([...top]);
  });

  it('returns one opaque RGBA color per scalar', () => {
    const colors = colorizeLogError(new Float32Array([-16, -8, 0]));
    expect(colors).toHaveLength(12);
    expect([colors[3], colors[7], colors[11]]).toEqual([1, 1, 1]);
  });

  it('rejects a non-finite field', () => {
    expect(() => colorizeLogError(new Float32Array([Number.NaN]))).toThrow(/not finite/i);
  });
});
