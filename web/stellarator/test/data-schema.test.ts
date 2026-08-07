import { describe, expect, it } from 'vitest';
import {
  SUPPORTED_ORDERS,
  validateDiscretization,
  validateOrder,
  validateRestol,
  validateSolverDataset,
  validateSurface,
} from '../src/data-schema';

function validDataset() {
  return {
    schema: 1 as const,
    nsrc: 2,
    nrender: 3,
    nfaces: 1,
    grfError: 1e-3,
    elapsedMs: 12,
    positions: new Float32Array(9),
    logError: new Float32Array([-4, -3, -2]),
    triangles: new Uint32Array([0, 1, 2]),
  };
}

describe('validateSolverDataset', () => {
  it('accepts a finite, consistently shaped mesh', () => {
    expect(validateSolverDataset(validDataset())).toEqual(validDataset());
  });

  it('rejects a triangle index outside the render lattice', () => {
    const data = validDataset();
    data.triangles[2] = 3;
    expect(() => validateSolverDataset(data)).toThrow(/triangle index/i);
  });

  it('rejects non-finite field values', () => {
    const data = validDataset();
    data.logError[1] = Number.NaN;
    expect(() => validateSolverDataset(data)).toThrow(/logError/i);
  });
});

describe('validateDiscretization', () => {
  it('passes positive safe integer mp and np through unchanged', () => {
    expect(validateDiscretization(8, 24)).toEqual({ mp: 8, np: 24 });
  });

  it.each([[0, 24], [8, 0], [1.5, 24], [8, Number.NaN]])(
    'rejects invalid dimensions %s × %s',
    (mp, np) => {
      expect(() => validateDiscretization(mp, np)).toThrow(/mp and np/i);
    },
  );
});

describe('validateOrder', () => {
  const expectedOrders = [4, 6, 8, 10, 12, 14, 16] as const;

  it('exposes the even stellarator orders from 4 through 16', () => {
    expect(SUPPORTED_ORDERS).toEqual(expectedOrders);
  });

  it.each(expectedOrders)('accepts supported order %i', (order) => {
    expect(validateOrder(order)).toBe(order);
  });

  it.each([0, 2, 5, 18, Number.NaN])('rejects unsupported order %s', (order) => {
    expect(() => validateOrder(order)).toThrow(/order/i);
  });
});

describe('validateSurface', () => {
  it('accepts the two supported surfaces', () => {
    expect(validateSurface('builtin')).toBe('builtin');
    expect(validateSurface('w7x')).toBe('w7x');
  });

  it.each([['torus'], [''], ['W7X']])('rejects surface %s', (surface) => {
    expect(() => validateSurface(surface)).toThrow(/surface/i);
  });
});

describe('validateRestol', () => {
  it('passes a positive finite tolerance through', () => {
    expect(validateRestol(0.1)).toBe(0.1);
  });

  it.each([[0], [-0.1], [Number.NaN], [Number.POSITIVE_INFINITY]])(
    'rejects restol %s',
    (restol) => {
      expect(() => validateRestol(restol)).toThrow(/restol/i);
    },
  );
});
