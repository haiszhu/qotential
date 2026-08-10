import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath, pathToFileURL } from 'node:url';

const fixturePath = process.argv[2];
if (!fixturePath) {
  throw new Error('usage: wasm64-order12-parity.mjs <native-fixture.bin>');
}

const here = path.dirname(fileURLToPath(import.meta.url));
const webRoot = path.resolve(here, '..');
const moduleUrl = pathToFileURL(path.join(webRoot, 'public/wasm/solver64.js')).href;
const wasmPath = path.join(webRoot, 'public/wasm/solver64.wasm');
const { default: createSolver } = await import(moduleUrl);
const wasm = await createSolver({ wasmBinary: fs.readFileSync(wasmPath) });

function require(condition, message) {
  if (!condition) throw new Error(message);
}

function parseNativeFixture(filename) {
  const blob = fs.readFileSync(filename);
  require(blob.subarray(0, 8).toString('ascii') === 'STGRF001', 'bad fixture magic');
  const ntri = Number(blob.readBigInt64LE(8));
  const nsrc = Number(blob.readBigInt64LE(16));
  const nrender = Number(blob.readBigInt64LE(24));
  const nfaces = Number(blob.readBigInt64LE(32));
  const grfError = blob.readDoubleLE(40);
  require(ntri > 0 && nsrc > 0 && nrender > 0 && nfaces > 0,
    'native fixture has invalid dimensions');
  require(Number.isFinite(grfError), 'native fixture GRF is non-finite');

  const ubnOffset = 48 + 8 * 8 * nsrc;
  let ubnScale = 0;
  for (let i = 0; i < nsrc; ++i) {
    const value = blob.readDoubleLE(ubnOffset + 8 * i);
    require(Number.isFinite(value), `native fixture ubn[${i}] is non-finite`);
    ubnScale = Math.max(ubnScale, Math.abs(value));
  }
  require(ubnScale > 0, 'native fixture ubn scale is not positive');

  const uOffset = 48 + 8 * 9 * nsrc;
  const u = new Float64Array(nsrc);
  for (let i = 0; i < nsrc; ++i) u[i] = blob.readDoubleLE(uOffset + 8 * i);
  require(u.every(Number.isFinite), 'native fixture u contains non-finite values');

  const triangleBytes = 8 * 3 * nfaces;
  const triangleOffset = blob.length - triangleBytes;
  require(triangleOffset >= uOffset + 8 * nsrc, 'native fixture is truncated');
  const triangles = new BigInt64Array(3 * nfaces);
  for (let i = 0; i < triangles.length; ++i) {
    triangles[i] = blob.readBigInt64LE(triangleOffset + 8 * i);
  }
  return { ntri, nsrc, nrender, nfaces, grfError, ubnScale, u, triangles };
}

function pointerArgument(offset) {
  require(Number.isSafeInteger(offset) && offset > 0, `invalid heap offset: ${offset}`);
  return BigInt(offset);
}

function allocateF64(count, name) {
  const offset = wasm._malloc(8 * count);
  require(offset !== 0, `${name}: WASM allocation failed`);
  return { offset, pointer: pointerArgument(offset) };
}

function createSimplexPrecomp(order) {
  const nquad = order + 2;
  const korder = order - 1;
  const kpols = order * (order + 1) / 2;
  const counts = [nquad, nquad, nquad * nquad, nquad,
    nquad * nquad, kpols * kpols, kpols * kpols];
  const buffers = counts.map((count, index) => allocateF64(count, `simplex-${index}`));
  const status = wasm._solver_simplex_precomp(
    BigInt(nquad), BigInt(korder), BigInt(kpols),
    ...buffers.map((buffer) => buffer.pointer),
  );
  require(status === 0,
    `simplex precomputation failed: ${status}/${wasm._solver_last_error()}`);
  return {
    pointers: buffers.map((buffer) => buffer.pointer),
    dispose() {
      for (let i = buffers.length - 1; i >= 0; --i) wasm._free(buffers[i].pointer);
    },
  };
}

function copyF64(name, count) {
  const buffer = allocateF64(count, name);
  try {
    const status = wasm[name](buffer.pointer, BigInt(count));
    require(status === 0, `${name} failed: ${status}/${wasm._solver_last_error()}`);
    return new Float64Array(wasm.HEAPU8.buffer, buffer.offset, count).slice();
  } finally {
    wasm._free(buffer.pointer);
  }
}

function copyI64(name, count) {
  const buffer = allocateF64(count, name);
  try {
    const status = wasm[name](buffer.pointer, BigInt(count));
    require(status === 0, `${name} failed: ${status}/${wasm._solver_last_error()}`);
    return new BigInt64Array(wasm.HEAPU8.buffer, buffer.offset, count).slice();
  } finally {
    wasm._free(buffer.pointer);
  }
}

const native = parseNativeFixture(fixturePath);
const refs = createSimplexPrecomp(12);
const events = [];
wasm.onSolverProgress = (stage, current, total, aux0, aux1, value) => {
  events.push({
    stage: Number(stage), current: Number(current), total: Number(total),
    aux0: Number(aux0), aux1: Number(aux1), value,
  });
  if (stage === 1 || stage === 2 || stage === 3 || stage === 4 ||
      stage === 9 || stage === 10 || stage === 11 ||
      ((stage === 6 || stage === 8) && current === total)) {
    console.log('WASM64_ORDER12_PROGRESS', Number(stage),
      `${String(current)}/${String(total)}`, value);
  }
};

const started = performance.now();
let status;
try {
  status = wasm._solver_run(
    48n, 144n, 12n, 0n, 0.1, 1n, 1e-12,
    ...refs.pointers,
  );
} finally {
  delete wasm.onSolverProgress;
}
require(status === 0, `order-12 solve failed: ${status}/${wasm._solver_last_error()}`);
const elapsedSeconds = (performance.now() - started) / 1000;

const nsrc = Number(wasm._solver_result_nsrc());
const nrender = Number(wasm._solver_result_nrender());
const nfaces = Number(wasm._solver_result_ntriangles());
const grfError = wasm._solver_result_grf_error();
require(nsrc === native.nsrc && nrender === native.nrender && nfaces === native.nfaces,
  `topology dimensions differ: wasm=${nsrc}/${nrender}/${nfaces} ` +
  `native=${native.nsrc}/${native.nrender}/${native.nfaces}`);

const u = copyF64('_solver_copy_u', nsrc);
const triangles = copyI64('_solver_copy_render_triangles', 3 * nfaces);
require(triangles.every((value, index) => value === native.triangles[index]),
  'render topology differs from native');

let maxDifference = 0;
for (let i = 0; i < nsrc; ++i) {
  require(Number.isFinite(u[i]), `WASM u[${i}] is non-finite`);
  maxDifference = Math.max(maxDifference, Math.abs(u[i] - native.u[i]));
}
// Use the same physical scale as the GRF definition.  Normalizing the
// cancellation residual by itself makes tiny sign changes look artificially
// large and does not measure the boundary-data-relative field error.
const normalizedFieldDifference = maxDifference / native.ubnScale;
const heapBytes = wasm.HEAPU8.buffer.byteLength;
require(normalizedFieldDifference <= 1e-9,
  `normalized u difference ${normalizedFieldDifference} exceeds 1e-9`);
require(Number.isFinite(grfError) && grfError < 2e-9,
  `WASM GRF ${grfError} is outside the order-12 envelope`);
require(Math.abs(grfError - native.grfError) < 1.2e-9,
  `WASM/native GRF difference ${Math.abs(grfError-native.grfError)} exceeds 1.2e-9`);
require(Math.abs(native.grfError - 8.5129424870771411e-10) < 1e-15,
  `native GRF anchor changed: ${native.grfError}`);
require(heapBytes < 17179869184,
  `WASM heap ${heapBytes} reached the 16 GiB limit`);

for (const stage of [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]) {
  require(events.some((event) => event.stage === stage), `missing progress stage ${stage}`);
}
for (const stage of [4, 6, 8]) {
  const complete = events.filter((event) =>
    event.stage === stage && event.current === event.total);
  require(complete.length === 1, `stage ${stage} did not complete exactly once`);
}

console.log('WASM64_ORDER12_PARITY_OK',
  `ntri=${native.ntri}`, `nsrc=${nsrc}`, `nrender=${nrender}`,
  `normalizedUByUbn=${normalizedFieldDifference}`,
  `nativeGrf=${native.grfError}`, `wasmGrf=${grfError}`,
  `heap=${heapBytes}`, `seconds=${elapsedSeconds}`);

wasm._solver_clear();
refs.dispose();
