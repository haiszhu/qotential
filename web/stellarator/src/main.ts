import './style.css';
import { runAdvisory } from './run-advisory';
import {
  validateDiscretization,
  validateFmmTolerance,
  validateKernel,
  validateOrder,
  validateRestol,
  validateSolverDataset,
  validateSurface,
  type WorkerResponse,
} from './data-schema';
import { prepareGpuMesh, StellaratorRenderer } from './renderer';
import { createLogBuffer, shouldFollow } from './log-view';

function element<T extends HTMLElement>(selector: string): T {
  const value = document.querySelector<T>(selector);
  if (!value) throw new Error(`Missing required element: ${selector}`);
  return value;
}

const canvas = element<HTMLCanvasElement>('#stellarator-canvas');
const runButton = element<HTMLButtonElement>('#run-solver');
const status = element<HTMLElement>('#run-status');
const statusDot = element<HTMLElement>('#status-dot');
const grfValue = element<HTMLElement>('#grf-value');
const sourceValue = element<HTMLElement>('#source-value');
const faceValue = element<HTMLElement>('#face-value');
const timeValue = element<HTMLElement>('#time-value');
const placeholder = element<HTMLElement>('#canvas-placeholder');
const mpInput = element<HTMLInputElement>('#mp-input');
const npInput = element<HTMLInputElement>('#np-input');
const orderInput = element<HTMLSelectElement>('#order-input');
const surfaceInput = element<HTMLSelectElement>('#surface-input');
const restolInput = element<HTMLInputElement>('#restol-input');
const kernelInput = element<HTMLSelectElement>('#kernel-input');
const fmmToleranceInput = element<HTMLSelectElement>('#fmm-tolerance-input');
const advisory = element<HTMLParagraphElement>('#run-advisory');
const solverLog = element<HTMLElement>('#solver-log');
// Authoritative line record; kept 1:1 with the rendered children so the
// clear/append semantics verified in log-view.test.ts are the page's own.
const logBuffer = createLogBuffer();

function clearLog(): void {
  logBuffer.clear();
  solverLog.replaceChildren();
}

// Append one child per line so role="log" stays incremental and assistive
// technology is not made to re-announce the whole transcript each event.
function appendLog(line: string): void {
  logBuffer.append(line);
  const follow = shouldFollow(
    solverLog.scrollTop, solverLog.clientHeight, solverLog.scrollHeight);
  const node = document.createElement('p');
  node.className = 'log-line';
  if (line.startsWith('[error]')) node.dataset.kind = 'error';
  else if (line.startsWith('[cancelled]')) node.dataset.kind = 'cancelled';
  node.textContent = line;
  solverLog.appendChild(node);
  if (follow) solverLog.scrollTop = solverLog.scrollHeight;
}

// Advice only: never change, disable, or reject the chosen discretization.
function refreshAdvisory(): void {
  restolInput.disabled = surfaceInput.value !== 'w7x';
  const text = runAdvisory(
    surfaceInput.value as 'builtin' | 'w7x',
    Number(orderInput.value),
    Number(mpInput.value),
    Number(npInput.value),
    Number(restolInput.value),
    validateKernel(kernelInput.value),
    validateFmmTolerance(Number(fmmToleranceInput.value)),
  );
  advisory.textContent = text;
  advisory.hidden = text === '';
}

for (const control of [surfaceInput, orderInput, mpInput, npInput, restolInput,
  kernelInput, fmmToleranceInput]) {
  control.addEventListener('change', refreshAdvisory);
  control.addEventListener('input', refreshAdvisory);
}
refreshAdvisory();

let renderer: StellaratorRenderer | undefined;
let requestId = 0;
let activeRequest = 0;
let running = false;
let worker = createWorker();

function setState(kind: 'idle' | 'running' | 'ready' | 'error', message: string): void {
  status.textContent = message;
  statusDot.dataset.state = kind;
  running = kind === 'running';
  runButton.disabled = !renderer;
  runButton.textContent = running ? 'Cancel solve' : 'Run Fortran solver';
}

async function initialize(): Promise<void> {
  try {
    renderer = await StellaratorRenderer.create(canvas);
    setState('idle', 'Ready — computation runs locally in your browser');
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    placeholder.textContent = `${message}. Use a current Chrome or Edge build with WebGPU enabled.`;
    setState('error', 'WebGPU unavailable');
  }
}

runButton.addEventListener('click', () => {
  if (!renderer) return;
  if (running) {
    worker.terminate();
    worker = createWorker();
    activeRequest = 0;
    placeholder.hidden = false;
    placeholder.textContent = 'The browser solve was cancelled. No partial result was rendered.';
    appendLog('[cancelled] solve terminated by user');
    setState('idle', 'Cancelled — ready to run again');
    return;
  }
  activeRequest = ++requestId;
  clearLog();
  let discretization: { mp: number; np: number };
  let order;
  let surface;
  let restol;
  let kernel;
  let fmmTolerance;
  try {
    discretization = validateDiscretization(Number(mpInput.value), Number(npInput.value));
    order = validateOrder(Number(orderInput.value));
    surface = validateSurface(surfaceInput.value);
    restol = surface === 'w7x' ? validateRestol(Number(restolInput.value)) : 0.1;
    kernel = validateKernel(kernelInput.value);
    fmmTolerance = validateFmmTolerance(Number(fmmToleranceInput.value));
  } catch (error) {
    const message = error instanceof Error ? error.message : String(error);
    placeholder.hidden = false;
    placeholder.textContent = message;
    appendLog(`[error] ${message}`);
    setState('error', 'Invalid discretization');
    return;
  }
  placeholder.hidden = true;
  setState('running', 'Starting WebAssembly worker');
  const wasm32ModuleUrl = new URL(
    `${import.meta.env.BASE_URL}wasm/solver.js`, document.baseURI).href;
  const wasm64ModuleUrl = new URL(
    `${import.meta.env.BASE_URL}wasm/solver64.js`, document.baseURI).href;
  worker.postMessage({
    type: 'run', requestId: activeRequest, wasm32ModuleUrl, wasm64ModuleUrl,
    ...discretization, order,
    surface, restol, kernel, fmmTolerance,
  });
});

function createWorker(): Worker {
  const next = new Worker(new URL('./solver-worker.ts', import.meta.url), { type: 'module' });
  next.onmessage = (event: MessageEvent<WorkerResponse>) => {
    const message = event.data;
    if (message.requestId !== activeRequest) return;
    if (message.type === 'progress') {
      setState('running', message.phase);
      return;
    }
    if (message.type === 'log') {
      appendLog(message.line);
      return;
    }
    if (message.type === 'error') {
      placeholder.hidden = false;
      placeholder.textContent = message.message;
      appendLog(`[error] ${message.message}`);
      next.terminate();
      worker = createWorker();
      activeRequest = 0;
      setState('error', 'Solver failed');
      return;
    }

    try {
      const data = validateSolverDataset(message.data);
      renderer!.setMesh(prepareGpuMesh(data.positions, data.logError, data.triangles));
      grfValue.textContent = data.grfError.toExponential(6);
      sourceValue.textContent = data.nsrc.toLocaleString();
      faceValue.textContent = data.nfaces.toLocaleString();
      timeValue.textContent = `${(data.elapsedMs / 1000).toFixed(2)} s`;
      setState('ready', 'Complete — drag to rotate, scroll to zoom');
    } catch (error) {
      // The solve may have succeeded and already logged its GRF/completed lines,
      // but result validation, mesh preparation, or WebGPU upload can still fail
      // here.  Append the error so the log does not end on a false success.
      const message = error instanceof Error ? error.message : String(error);
      placeholder.hidden = false;
      placeholder.textContent = message;
      appendLog(`[error] ${message}`);
      setState('error', 'Invalid solver result');
    }
  };
  next.onerror = (event) => {
    const message = event.message || 'The WebAssembly worker stopped unexpectedly.';
    placeholder.hidden = false;
    placeholder.textContent = message;
    appendLog(`[error] ${message}`);
    next.terminate();
    worker = createWorker();
    activeRequest = 0;
    setState('error', 'Worker error');
  };
  return next;
}

void initialize();
