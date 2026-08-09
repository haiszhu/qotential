export interface SolverWasmModule {
  HEAPU8: Uint8Array;
  _malloc(bytes: number): number;
  _free(pointer: number): void;
  _solver_clear(): void;
  _solver_simplex_precomp(
    nquad: bigint, korder: bigint, kpols: bigint,
    tgl: number, wgl: number, Dgl: number, wBclag: number,
    Legmat: number, umatr: number, vmatr: number,
  ): number;
  _solver_run(
    mp: bigint, np: bigint, order: bigint, surface: bigint, restol: number,
    kernel: bigint, fmmTolerance: number,
    tgl: number, wgl: number, Dgl: number, wBclag: number,
    Legmat: number, umatr: number, vmatr: number,
  ): number;
  _solver_result_nsrc(): bigint | number;
  _solver_result_nrender(): bigint | number;
  _solver_result_ntriangles(): bigint | number;
  _solver_result_grf_error(): number;
  _solver_last_error(): number;
  _solver_copy_render_xyz(pointer: number, capacity: bigint): number;
  _solver_copy_render_log_error(pointer: number, capacity: bigint): number;
  _solver_copy_render_triangles(pointer: number, capacity: bigint): number;
  // Installed per solve by the Worker; invoked from the EM_JS progress bridge.
  // With -sWASM_BIGINT=1 the four i64 payloads arrive as bigint.
  onSolverProgress?: (
    stage: number,
    current: bigint,
    total: bigint,
    aux0: bigint,
    aux1: bigint,
    value: number,
  ) => void;
}

export interface SimplexPrecomp {
  readonly tgl: number;
  readonly wgl: number;
  readonly Dgl: number;
  readonly wBclag: number;
  readonly Legmat: number;
  readonly umatr: number;
  readonly vmatr: number;
  dispose(): void;
}

export interface RawSolverResult {
  nsrc: number;
  nrender: number;
  nfaces: number;
  grfError: number;
  positions: Float64Array;
  logError: Float64Array;
  triangles: BigInt64Array;
}

export class SolverStatusError extends Error {
  constructor(
    readonly phase: string,
    readonly status: number,
    readonly solverError: number,
  ) {
    super(`${phase} failed: status=${status}, solver error=${solverError}`);
    this.name = 'SolverStatusError';
  }
}

export function requireStatus(
  phase: string,
  status: number,
  lastError: () => number,
): void {
  if (status !== 0) {
    throw new SolverStatusError(phase, status, lastError());
  }
}

export function createSimplexPrecomp(
  wasm: SolverWasmModule,
  order: number,
): SimplexPrecomp {
  const nquad = order + 2;
  const korder = order - 1;
  const kpols = order * (order + 1) / 2;
  const pointers: number[] = [];
  const allocate = (count: number, name: string): number => {
    const pointer = wasm._malloc(8 * count);
    if (!pointer) throw new Error(`${name}: WASM allocation failed`);
    pointers.push(pointer);
    return pointer;
  };

  try {
    const tgl = allocate(nquad, 'tgl');
    const wgl = allocate(nquad, 'wgl');
    const Dgl = allocate(nquad * nquad, 'Dgl');
    const wBclag = allocate(nquad, 'w_bclag');
    const Legmat = allocate(nquad * nquad, 'Legmat');
    const umatr = allocate(kpols * kpols, 'umatr');
    const vmatr = allocate(kpols * kpols, 'vmatr');
    requireStatus('simplex precomputation', wasm._solver_simplex_precomp(
      BigInt(nquad), BigInt(korder), BigInt(kpols),
      tgl, wgl, Dgl, wBclag, Legmat, umatr, vmatr,
    ), () => wasm._solver_last_error());

    let disposed = false;
    return {
      tgl, wgl, Dgl, wBclag, Legmat, umatr, vmatr,
      dispose() {
        if (disposed) return;
        disposed = true;
        for (let i = pointers.length - 1; i >= 0; --i) wasm._free(pointers[i]!);
      },
    };
  } catch (error) {
    for (let i = pointers.length - 1; i >= 0; --i) wasm._free(pointers[i]!);
    throw error;
  }
}

export async function loadSolverModule(moduleUrl: string): Promise<SolverWasmModule> {
  const imported = await import(/* @vite-ignore */ moduleUrl) as {
    default: (options?: Record<string, unknown>) => Promise<SolverWasmModule>;
  };
  return imported.default();
}

function copyF64(
  wasm: SolverWasmModule,
  count: number,
  copy: (pointer: number, capacity: bigint) => number,
  name: string,
): Float64Array {
  const pointer = wasm._malloc(8 * count);
  if (!pointer) throw new Error(`${name}: WASM allocation failed`);
  try {
    requireStatus(name, copy(pointer, BigInt(count)), () => wasm._solver_last_error());
    return new Float64Array(wasm.HEAPU8.buffer, pointer, count).slice();
  } finally {
    wasm._free(pointer);
  }
}

function copyI64(
  wasm: SolverWasmModule,
  count: number,
  copy: (pointer: number, capacity: bigint) => number,
  name: string,
): BigInt64Array {
  const pointer = wasm._malloc(8 * count);
  if (!pointer) throw new Error(`${name}: WASM allocation failed`);
  try {
    requireStatus(name, copy(pointer, BigInt(count)), () => wasm._solver_last_error());
    return new BigInt64Array(wasm.HEAPU8.buffer, pointer, count).slice();
  } finally {
    wasm._free(pointer);
  }
}

export function collectSolverResult(wasm: SolverWasmModule): RawSolverResult {
  const nsrc = Number(wasm._solver_result_nsrc());
  const nrender = Number(wasm._solver_result_nrender());
  const nfaces = Number(wasm._solver_result_ntriangles());
  if (nsrc <= 0 || nrender <= 0 || nfaces <= 0) {
    throw new Error(`solver returned invalid sizes: ${nsrc}/${nrender}/${nfaces}`);
  }
  return {
    nsrc,
    nrender,
    nfaces,
    grfError: wasm._solver_result_grf_error(),
    positions: copyF64(wasm, 3 * nrender,
      wasm._solver_copy_render_xyz.bind(wasm), 'copy positions'),
    logError: copyF64(wasm, nrender,
      wasm._solver_copy_render_log_error.bind(wasm), 'copy log error'),
    triangles: copyI64(wasm, 3 * nfaces,
      wasm._solver_copy_render_triangles.bind(wasm), 'copy triangles'),
  };
}
