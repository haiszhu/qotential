/// <reference lib="webworker" />

import type { WorkerRequest, WorkerResponse } from './data-schema';
import {
  collectSolverResult,
  createSimplexPrecomp,
  loadSolverModule,
  requireStatus,
  SolverStatusError,
  type SolverWasmModule,
} from './wasm-loader';
import { createProgressForwarder, withSolverProgress } from './progress-log';
import { createSolverRuntimeSelector } from './memory64';
import { toSolverDataset } from './worker-utils';

const scope: DedicatedWorkerGlobalScope = self as unknown as DedicatedWorkerGlobalScope;
let running = false;
let modulePromise: Promise<SolverWasmModule> | undefined;
let loadedUrl = '';
const selectRuntime = createSolverRuntimeSelector();

function send(message: WorkerResponse, transfer: Transferable[] = []): void {
  scope.postMessage(message, transfer);
}

function progress(requestId: number, phase: string): void {
  send({ type: 'progress', requestId, phase });
}

async function solverFor(url: string): Promise<SolverWasmModule> {
  if (!modulePromise || loadedUrl !== url) {
    loadedUrl = url;
    modulePromise = loadSolverModule(url);
  }
  return modulePromise;
}

scope.onmessage = async (event: MessageEvent<WorkerRequest>) => {
  const request = event.data;
  if (request.type !== 'run') return;
  if (running) {
    send({ type: 'error', requestId: request.requestId, message: 'The solver is already running.' });
    return;
  }

  running = true;
  const started = performance.now();
  let wasm: SolverWasmModule | undefined;
  try {
    progress(request.requestId, 'Loading Fortran WebAssembly runtime');
    const runtime = await selectRuntime({
      wasm32ModuleUrl: request.wasm32ModuleUrl,
      wasm64ModuleUrl: request.wasm64ModuleUrl,
    });
    send({ type: 'log', requestId: request.requestId, line: runtime.logLine });
    wasm = await solverFor(runtime.moduleUrl);
    const lastError = () => wasm!._solver_last_error();

    // The generated core exposes one atomic solve. Keep the
    // UI honest about that boundary instead of presenting adapter state changes
    // as separately timed numerical phases.
    const kernelDescription = request.kernel === 'fmm'
      ? `FMM3D (eps=${request.fmmTolerance})`
      : 'direct kernel';
    progress(request.requestId,
      `Running ${request.surface} surface at order ${request.order}, ` +
      `${request.mp} × ${request.np}, ${kernelDescription}, and close correction`);
    // Real Fortran stages append log lines through this forwarder while the
    // atomic solve runs.  withSolverProgress installs the callback only for the
    // duration of _solver_run and clears it before any render buffers are read.
    const forwarder = createProgressForwarder({
      sink: (line) => send({ type: 'log', requestId: request.requestId, line }),
    });
    const solveModule = wasm;
    const refs = createSimplexPrecomp(solveModule, request.order, runtime.pointerBits);
    try {
      withSolverProgress(solveModule, forwarder, () => {
        requireStatus(`order-${request.order} ${request.surface} solve`,
          solveModule._solver_run(
            BigInt(request.mp), BigInt(request.np), BigInt(request.order),
            BigInt(request.surface === 'w7x' ? 1 : 0), request.restol,
            BigInt(request.kernel === 'fmm' ? 1 : 0), request.fmmTolerance,
            refs.tgl, refs.wgl, refs.Dgl, refs.wBclag,
            refs.Legmat, refs.umatr, refs.vmatr,
          ), lastError);
      });
    } finally {
      refs.dispose();
    }
    progress(request.requestId, 'Preparing the render mesh');
    // Measure elapsed AFTER result assembly, matching the pre-progress metric
    // scope and the design's "after _solver_run and collectSolverResult() both
    // succeed" — the buffer copies count toward the reported browser-solve time.
    const rawResult = collectSolverResult(wasm, runtime.pointerBits);
    const elapsedMs = performance.now() - started;
    const data = toSolverDataset(rawResult, elapsedMs);
    // One extra truthful line after both the solve and result assembly succeed;
    // the GRF line itself arrives from the Fortran stage-11 callback above.
    send({ type: 'log', requestId: request.requestId,
      line: `[result] completed in ${(elapsedMs / 1000).toFixed(1)} s` });
    send({ type: 'result', requestId: request.requestId, data }, [
      data.positions.buffer,
      data.logError.buffer,
      data.triangles.buffer,
    ]);
  } catch (error) {
    const message = error instanceof SolverStatusError
      ? `status=${error.status}, solver error=${error.solverError}`
      : error instanceof Error ? error.message : String(error);
    send({ type: 'error', requestId: request.requestId, message });
  } finally {
    wasm?._solver_clear();
    running = false;
  }
};

export {};
