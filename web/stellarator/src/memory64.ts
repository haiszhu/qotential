import type { PointerBits } from './wasm-loader';

const memory64Probe = Uint8Array.from([
  0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00,
  0x05, 0x07, 0x01, 0x05, 0x80, 0x10, 0x80, 0x80, 0x10,
]);

export interface Memory64ProbeRuntime {
  validate(bytes: BufferSource): boolean;
  instantiate(bytes: BufferSource): Promise<unknown>;
}

export interface SolverModuleUrls {
  wasm32ModuleUrl: string;
  wasm64ModuleUrl: string;
}

export interface SolverRuntimeSelection {
  moduleUrl: string;
  pointerBits: PointerBits;
  logLine: string;
}

export async function detectMemory64(
  runtime: Memory64ProbeRuntime = WebAssembly,
): Promise<boolean> {
  if (!runtime.validate(memory64Probe)) return false;
  await runtime.instantiate(memory64Probe);
  return true;
}

export async function selectSolverRuntime(
  urls: SolverModuleUrls,
  probe: () => Promise<boolean> = detectMemory64,
): Promise<SolverRuntimeSelection> {
  if (await probe()) {
    return {
      moduleUrl: urls.wasm64ModuleUrl,
      pointerBits: 64,
      logLine: '[runtime] WebAssembly Memory64 enabled',
    };
  }
  return {
    moduleUrl: urls.wasm32ModuleUrl,
    pointerBits: 32,
    logLine: '[runtime] Memory64 unavailable; using wasm32',
  };
}

export function createSolverRuntimeSelector(
  probe: () => Promise<boolean> = detectMemory64,
): (urls: SolverModuleUrls) => Promise<SolverRuntimeSelection> {
  let selected: Promise<SolverRuntimeSelection> | undefined;
  return (urls) => {
    selected ??= selectSolverRuntime(urls, probe);
    return selected;
  };
}
