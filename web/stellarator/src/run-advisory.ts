import type { SolverOrder, SolverSurface } from './data-schema';

/**
 * The one browser configuration whose W7-X source count has actually been
 * measured and frozen as a native parity fixture. Everything else is
 * data-dependent: the curvature criterion decides the chart count during real
 * geometry construction, so N cannot be predicted from `restol` beforehand.
 */
export const MEASURED_W7X_REFERENCE = {
  mp: 12,
  np: 36,
  order: 6,
  restol: 0.1,
  nsrc: 115_500,
} as const;

const NODES = MEASURED_W7X_REFERENCE.nsrc.toLocaleString('en-US');

const REFERENCE_PHRASE =
  `order ${MEASURED_W7X_REFERENCE.order} at ` +
  `${MEASURED_W7X_REFERENCE.mp} × ${MEASURED_W7X_REFERENCE.np}, ` +
  `restol ${MEASURED_W7X_REFERENCE.restol} (N=${NODES}) is the ` +
  `measured reference`;

function isMeasuredReference(
  order: number, mp: number, np: number, restol: number,
): boolean {
  return order === MEASURED_W7X_REFERENCE.order &&
    mp === MEASURED_W7X_REFERENCE.mp &&
    np === MEASURED_W7X_REFERENCE.np &&
    restol === MEASURED_W7X_REFERENCE.restol;
}

/**
 * Advice only. This never changes, disables, or rejects the chosen values —
 * the discretization is the user's decision.
 */
export function runAdvisory(
  surface: SolverSurface,
  order: SolverOrder | number,
  mp: number,
  np: number,
  restol: number,
): string {
  if (surface !== 'w7x') return '';

  if (order === 4) {
    return `W7-X at order 4 has no frozen W7-X parity reference; ` +
      `${REFERENCE_PHRASE}. The direct sum is O(N²), so solve time grows ` +
      `quadratically with the source count.`;
  }

  if (isMeasuredReference(order, mp, np, restol)) {
    return `This is the measured W7-X reference (N=${NODES}): ` +
      `${NODES} source nodes. The direct sum is O(N²).`;
  }

  return `W7-X source node count is data-dependent: curvature refinement ` +
    `decides the chart count, so N is only known after real geometry is ` +
    `built. For scale, ${REFERENCE_PHRASE}.`;
}
