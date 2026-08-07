import { validateSolverDataset, type SolverDataset } from './data-schema';
import type { RawSolverResult } from './wasm-loader';

export function toSolverDataset(raw: RawSolverResult, elapsedMs: number): SolverDataset {
  const triangles = new Uint32Array(raw.triangles.length);
  for (let i = 0; i < raw.triangles.length; ++i) {
    // The public WASM adapter deliberately exports WebGPU-ready zero-based indices.
    const index = Number(raw.triangles[i]);
    if (!Number.isSafeInteger(index) || index < 0 || index >= raw.nrender) {
      throw new Error(`triangle index ${raw.triangles[i]} at ${i} is invalid`);
    }
    triangles[i] = index;
  }
  return validateSolverDataset({
    schema: 1,
    nsrc: raw.nsrc,
    nrender: raw.nrender,
    nfaces: raw.nfaces,
    grfError: raw.grfError,
    elapsedMs,
    positions: Float32Array.from(raw.positions),
    logError: Float32Array.from(raw.logError),
    triangles,
  });
}
