export interface SolverDataset {
  schema: 1;
  nsrc: number;
  nrender: number;
  nfaces: number;
  grfError: number;
  elapsedMs: number;
  positions: Float32Array;
  logError: Float32Array;
  triangles: Uint32Array;
}

export const SUPPORTED_ORDERS = [4, 6, 8, 10, 12, 14, 16] as const;
export type SolverOrder = (typeof SUPPORTED_ORDERS)[number];

export const SUPPORTED_SURFACES = ['builtin', 'w7x'] as const;
export type SolverSurface = (typeof SUPPORTED_SURFACES)[number];

export const SUPPORTED_KERNELS = ['fmm', 'direct'] as const;
export type SolverKernel = (typeof SUPPORTED_KERNELS)[number];

export const SUPPORTED_FMM_TOLERANCES =
  [1e-3, 1e-6, 1e-9, 1e-12, 1e-15] as const;
export type FmmTolerance = (typeof SUPPORTED_FMM_TOLERANCES)[number];

export type WorkerRequest = {
  type: 'run';
  requestId: number;
  moduleUrl: string;
  mp: number;
  np: number;
  order: SolverOrder;
  surface: SolverSurface;
  restol: number;
  kernel: SolverKernel;
  fmmTolerance: FmmTolerance;
};

export type WorkerResponse =
  | { type: 'progress'; requestId: number; phase: string }
  | { type: 'log'; requestId: number; line: string }
  | { type: 'result'; requestId: number; data: SolverDataset }
  | { type: 'error'; requestId: number; message: string };

function requireCondition(condition: unknown, message: string): asserts condition {
  if (!condition) throw new Error(message);
}

function requireFinite(name: string, values: Float32Array): void {
  for (let i = 0; i < values.length; ++i) {
    requireCondition(Number.isFinite(values[i]), `${name}[${i}] is not finite`);
  }
}

export function validateDiscretization(mp: number, np: number): { mp: number; np: number } {
  requireCondition(
    Number.isSafeInteger(mp) && mp > 0 && Number.isSafeInteger(np) && np > 0,
    'mp and np must be positive integers',
  );
  return { mp, np };
}

export function validateOrder(order: number): SolverOrder {
  requireCondition(
    SUPPORTED_ORDERS.includes(order as SolverOrder),
    `order must be one of ${SUPPORTED_ORDERS.join(', ')}`,
  );
  return order as SolverOrder;
}

export function validateSurface(surface: string): SolverSurface {
  requireCondition(
    SUPPORTED_SURFACES.includes(surface as SolverSurface),
    `surface must be one of ${SUPPORTED_SURFACES.join(', ')}`,
  );
  return surface as SolverSurface;
}

export function validateRestol(restol: number): number {
  requireCondition(
    Number.isFinite(restol) && restol > 0,
    'restol must be a positive finite number',
  );
  return restol;
}

export function validateKernel(kernel: string): SolverKernel {
  requireCondition(
    SUPPORTED_KERNELS.includes(kernel as SolverKernel),
    `kernel must be one of ${SUPPORTED_KERNELS.join(', ')}`,
  );
  return kernel as SolverKernel;
}

export function validateFmmTolerance(tolerance: number): FmmTolerance {
  requireCondition(
    SUPPORTED_FMM_TOLERANCES.includes(tolerance as FmmTolerance),
    `FMM tolerance must be one of ${SUPPORTED_FMM_TOLERANCES.join(', ')}`,
  );
  return tolerance as FmmTolerance;
}

export function validateSolverDataset(data: SolverDataset): SolverDataset {
  requireCondition(data.schema === 1, 'unsupported solver dataset schema');
  requireCondition(Number.isInteger(data.nsrc) && data.nsrc > 0, 'invalid source count');
  requireCondition(Number.isInteger(data.nrender) && data.nrender > 0, 'invalid render count');
  requireCondition(Number.isInteger(data.nfaces) && data.nfaces > 0, 'invalid face count');
  requireCondition(Number.isFinite(data.grfError) && data.grfError >= 0, 'invalid GRF error');
  requireCondition(Number.isFinite(data.elapsedMs) && data.elapsedMs >= 0, 'invalid elapsed time');
  requireCondition(data.positions instanceof Float32Array, 'positions must be Float32Array');
  requireCondition(data.logError instanceof Float32Array, 'logError must be Float32Array');
  requireCondition(data.triangles instanceof Uint32Array, 'triangles must be Uint32Array');
  requireCondition(data.positions.length === 3 * data.nrender, 'position shape mismatch');
  requireCondition(data.logError.length === data.nrender, 'logError shape mismatch');
  requireCondition(data.triangles.length === 3 * data.nfaces, 'triangle shape mismatch');
  requireFinite('positions', data.positions);
  requireFinite('logError', data.logError);
  for (let i = 0; i < data.triangles.length; ++i) {
    requireCondition(data.triangles[i] < data.nrender,
      `triangle index ${data.triangles[i]} at ${i} is outside ${data.nrender}`);
  }
  return data;
}
