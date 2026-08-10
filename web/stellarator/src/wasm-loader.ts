export interface SolverWasmModule {
  HEAPU8: Uint8Array;
  _malloc(bytes: number): number;
  _free(pointer: WasmPointer): void;
  _solver_clear(): void;
  _solver_simplex_precomp(
    nquad: bigint, korder: bigint, kpols: bigint,
    tgl: WasmPointer, wgl: WasmPointer, Dgl: WasmPointer,
    wBclag: WasmPointer, Legmat: WasmPointer,
    umatr: WasmPointer, vmatr: WasmPointer,
  ): number;
  _solver_run(
    mp: bigint, np: bigint, order: bigint, surface: bigint, restol: number,
    kernel: bigint, fmmTolerance: number,
    tgl: WasmPointer, wgl: WasmPointer, Dgl: WasmPointer,
    wBclag: WasmPointer, Legmat: WasmPointer,
    umatr: WasmPointer, vmatr: WasmPointer,
  ): number;
  _solver_result_nsrc(): bigint | number;
  _solver_result_nrender(): bigint | number;
  _solver_result_ntriangles(): bigint | number;
  _solver_result_grf_error(): number;
  _solver_last_error(): number;
  _solver_copy_render_xyz(pointer: WasmPointer, capacity: bigint): number;
  _solver_copy_render_log_error(pointer: WasmPointer, capacity: bigint): number;
  _solver_copy_render_triangles(pointer: WasmPointer, capacity: bigint): number;
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

export type WasmPointer = number | bigint;
export type PointerBits = 32 | 64;

export function pointerArgument(offset: number, bits: PointerBits): WasmPointer {
  if (!Number.isSafeInteger(offset) || offset < 0) {
    throw new RangeError(`invalid WASM heap offset: ${offset}`);
  }
  return bits === 64 ? BigInt(offset) : offset;
}

export function checkedHeapOffset(pointer: WasmPointer): number {
  if (typeof pointer === 'number') {
    if (!Number.isSafeInteger(pointer) || pointer < 0) {
      throw new RangeError(`invalid WASM pointer result: ${pointer}`);
    }
    return pointer;
  }
  if (pointer < 0n || pointer > BigInt(Number.MAX_SAFE_INTEGER)) {
    throw new RangeError(`wasm64 pointer is outside the safe offset range: ${pointer}`);
  }
  return Number(pointer);
}

export function requireHeapRange(
  wasm: Pick<SolverWasmModule, 'HEAPU8'>,
  offset: number,
  byteLength: number,
): void {
  if (!Number.isSafeInteger(byteLength) || byteLength < 0 ||
      !Number.isSafeInteger(offset) || offset < 0 ||
      offset > wasm.HEAPU8.buffer.byteLength ||
      byteLength > wasm.HEAPU8.buffer.byteLength - offset) {
    throw new RangeError(
      `WASM heap range ${offset}+${byteLength} exceeds ${wasm.HEAPU8.buffer.byteLength}`,
    );
  }
}

export interface SimplexPrecomp {
  readonly tgl: WasmPointer;
  readonly wgl: WasmPointer;
  readonly Dgl: WasmPointer;
  readonly wBclag: WasmPointer;
  readonly Legmat: WasmPointer;
  readonly umatr: WasmPointer;
  readonly vmatr: WasmPointer;
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
  pointerBits: PointerBits = 32,
): SimplexPrecomp {
  const nquad = order + 2;
  const korder = order - 1;
  const kpols = order * (order + 1) / 2;
  const pointers: WasmPointer[] = [];
  const allocate = (count: number, name: string): WasmPointer => {
    const offset = wasm._malloc(8 * count);
    if (!offset) throw new Error(`${name}: WASM allocation failed`);
    const pointer = pointerArgument(offset, pointerBits);
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
  pointerBits: PointerBits,
  count: number,
  copy: (pointer: WasmPointer, capacity: bigint) => number,
  name: string,
): Float64Array {
  const offset = wasm._malloc(8 * count);
  if (!offset) throw new Error(`${name}: WASM allocation failed`);
  const pointer = pointerArgument(offset, pointerBits);
  try {
    requireStatus(name, copy(pointer, BigInt(count)), () => wasm._solver_last_error());
    requireHeapRange(wasm, offset, 8 * count);
    return new Float64Array(wasm.HEAPU8.buffer, offset, count).slice();
  } finally {
    wasm._free(pointer);
  }
}

function copyI64(
  wasm: SolverWasmModule,
  pointerBits: PointerBits,
  count: number,
  copy: (pointer: WasmPointer, capacity: bigint) => number,
  name: string,
): BigInt64Array {
  const offset = wasm._malloc(8 * count);
  if (!offset) throw new Error(`${name}: WASM allocation failed`);
  const pointer = pointerArgument(offset, pointerBits);
  try {
    requireStatus(name, copy(pointer, BigInt(count)), () => wasm._solver_last_error());
    requireHeapRange(wasm, offset, 8 * count);
    return new BigInt64Array(wasm.HEAPU8.buffer, offset, count).slice();
  } finally {
    wasm._free(pointer);
  }
}

export function collectSolverResult(
  wasm: SolverWasmModule,
  pointerBits: PointerBits = 32,
): RawSolverResult {
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
    positions: copyF64(wasm, pointerBits, 3 * nrender,
      wasm._solver_copy_render_xyz.bind(wasm), 'copy positions'),
    logError: copyF64(wasm, pointerBits, nrender,
      wasm._solver_copy_render_log_error.bind(wasm), 'copy log error'),
    triangles: copyI64(wasm, pointerBits, 3 * nfaces,
      wasm._solver_copy_render_triangles.bind(wasm), 'copy triangles'),
  };
}
