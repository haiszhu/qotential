import fs from 'node:fs';
import path from 'node:path';
import { createHash } from 'node:crypto';
import { fileURLToPath, pathToFileURL } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const webRoot = path.resolve(here, '..');
const artifact = process.env.WASM_ARTIFACT || 'solver';
const pointerBits = Number(process.env.WASM_POINTER_BITS || 32);
requirePointerBits(pointerBits);
const moduleUrl = pathToFileURL(path.join(webRoot, `public/wasm/${artifact}.js`)).href;
const wasmPath = path.join(webRoot, `public/wasm/${artifact}.wasm`);
const fixturePath = process.env.WASM_FIXTURE ||
  path.join(webRoot, 'fixtures/native/order4-direct.bin');
const order6FixturePath = path.join(webRoot, 'fixtures/native/order6-direct.bin');
const w7xFixturePath = path.join(webRoot, 'fixtures/native/w7x-order6-direct.bin');
const fmmFixtureDir = process.env.WASM_FMM_FIXTURE_DIR;
const requestedFmmCases = process.env.WASM_FMM_CASES
  ? new Set(process.env.WASM_FMM_CASES.split(',').filter(Boolean))
  : undefined;

if (!fs.statSync(wasmPath).size) throw new Error('empty solver.wasm');
const { default: createSolver } = await import(moduleUrl);
const wasm = await createSolver({ wasmBinary: fs.readFileSync(wasmPath) });

function requirePointerBits(bits) {
  if (bits !== 32 && bits !== 64) throw new Error(`WASM_POINTER_BITS must be 32 or 64: ${bits}`);
}

function pointerArgument(offset) {
  require(Number.isSafeInteger(offset) && offset >= 0, `invalid heap offset: ${offset}`);
  return pointerBits === 64 ? BigInt(offset) : offset;
}

function require(condition, message) {
  if (!condition) throw new Error(message);
}

function createSimplexPrecomp(order) {
  const nquad = order + 2;
  const korder = order - 1;
  const kpols = order * (order + 1) / 2;
  const counts = [nquad, nquad, nquad * nquad, nquad,
    nquad * nquad, kpols * kpols, kpols * kpols];
  const pointers = [];
  const offsets = [];
  try {
    for (const count of counts) {
      const pointer = wasm._malloc(8 * count);
      require(pointer !== 0, 'simplex precomputation malloc failed');
      offsets.push(pointer);
      pointers.push(pointerArgument(pointer));
    }
    const status = wasm._solver_simplex_precomp(
      BigInt(nquad), BigInt(korder), BigInt(kpols), ...pointers,
    );
    require(status === 0,
      `simplex precomputation failed: ${status}/${wasm._solver_last_error()}`);
    const values = offsets.map((offset, i) =>
      new Float64Array(wasm.HEAPU8.buffer, offset, counts[i]).slice());
    for (const field of values)
      require(field.every(Number.isFinite), 'non-finite simplex precomputation value');
    const tgl = values[0], wgl = values[1];
    require(tgl.every((value, i) => value > -1 && value < 1 &&
      (i === 0 || value > tgl[i - 1])), 'tgl is not strictly increasing in (-1,1)');
    require(Math.abs(wgl.reduce((sum, value) => sum + value, 0) - 2) < 1e-14,
      'Gauss weights do not sum to 2');
    let disposed = false;
    return {
      pointers, values,
      dispose() {
        if (disposed) return;
        disposed = true;
        for (let i = pointers.length - 1; i >= 0; --i) wasm._free(pointers[i]);
      },
    };
  } catch (error) {
    for (let i = pointers.length - 1; i >= 0; --i) wasm._free(pointers[i]);
    throw error;
  }
}

function solverRun(refs, mp, np, order, surface, restol,
                   kernel = 0n, fmmTolerance = 1e-3) {
  return wasm._solver_run(mp, np, order, surface, restol,
                          kernel, fmmTolerance, ...refs.pointers);
}

function parseFixture(filename) {
  const blob = fs.readFileSync(filename);
  require(blob.subarray(0, 8).toString('ascii') === 'STGRF001', 'bad fixture magic');
  const ntri = Number(blob.readBigInt64LE(8));
  const nsrc = Number(blob.readBigInt64LE(16));
  const nrender = Number(blob.readBigInt64LE(24));
  const nfaces = Number(blob.readBigInt64LE(32));
  const grfError = blob.readDoubleLE(40);
  let offset = 48;
  function f64(count) {
    const values = new Float64Array(count);
    for (let i = 0; i < count; ++i, offset += 8) values[i] = blob.readDoubleLE(offset);
    return values;
  }
  function i64(count) {
    const values = new BigInt64Array(count);
    for (let i = 0; i < count; ++i, offset += 8) values[i] = blob.readBigInt64LE(offset);
    return values;
  }
  const fields = {
    sx: f64(3 * nsrc), snx: f64(3 * nsrc), sw: f64(nsrc),
    ub: f64(nsrc), ubn: f64(nsrc), u: f64(nsrc),
    render_xyz: f64(3 * nrender), render_log_error: f64(nrender),
    render_triangles: i64(3 * nfaces),
  };
  require(offset === blob.length, 'fixture length mismatch');
  require(Number.isFinite(grfError), 'fixture GRF error is not finite');
  for (const [name, values] of Object.entries(fields)) {
    if (name === 'render_triangles') continue;
    require(values.every(Number.isFinite), `fixture ${name} contains non-finite values`);
  }
  return { ntri, nsrc, nrender, nfaces, grfError, fields };
}

const truth = parseFixture(fixturePath);
const order6Truth = parseFixture(order6FixturePath);
const w7xTruth = parseFixture(w7xFixturePath);
const refs4 = createSimplexPrecomp(4);
const refs4Repeat = createSimplexPrecomp(4);
const refs6 = createSimplexPrecomp(6);
for (let field = 0; field < refs4.values.length; ++field) {
  require(refs4.values[field].every(
    (value, index) => value === refs4Repeat.values[field][index]),
  `order-4 simplex precomputation field ${field} is not deterministic`);
}
refs4Repeat.dispose();

function callStage(name, ...args) {
  const status = wasm[name](...args);
  require(status === 0, `${name} failed: status=${status} last=${wasm._solver_last_error()}`);
}

require(solverRun(refs4, 0n, 36n, 4n, 0n, 0.1) === 105,
        'solver must reject a non-positive discretization');
wasm._solver_clear();
require(solverRun(refs4,
  BigInt(Number.MAX_SAFE_INTEGER), 52n, 4n, 0n, 0.1,
) === 105,
        'solver must reject overflowing derived result extents');
wasm._solver_clear();
require(solverRun(refs4, 4n, 12n, 5n, 0n, 0.1) === 106,
        'solver must reject an unsupported odd order');
wasm._solver_clear();

callStage('_solver_run', 12n, 36n, 6n, 0n, 0.1,
          0n, 1e-3, ...refs6.pointers);
const order6Nsrc = Number(wasm._solver_result_nsrc());
const order6Nrender = Number(wasm._solver_result_nrender());
const order6Nfaces = Number(wasm._solver_result_ntriangles());
require(order6Nsrc === order6Truth.nsrc,
        `order 6 source count=${order6Nsrc}, expected ${order6Truth.nsrc}`);
require(order6Nrender === order6Truth.nrender,
        `order 6 render count=${order6Nrender}, expected ${order6Truth.nrender}`);
require(order6Nfaces === order6Truth.nfaces,
        `order 6 face count=${order6Nfaces}, expected ${order6Truth.nfaces}`);
require(Math.abs(wasm._solver_result_grf_error() - order6Truth.grfError) < 1e-8,
        `order 6 GRF=${wasm._solver_result_grf_error()} differs from native`);
const order6RenderXyz = copyF64('_solver_copy_render_xyz', 3 * order6Nrender);
const order6RenderLog = copyF64('_solver_copy_render_log_error', order6Nrender);
const order6Triangles = copyI64('_solver_copy_render_triangles', 3 * order6Nfaces);
const order6DirectWasm = {
  grfError: wasm._solver_result_grf_error(),
  u: copyF64('_solver_copy_u', order6Nsrc),
  renderTriangles: order6Triangles,
};
require(order6Triangles.every(
  (value, index) => value === order6Truth.fields.render_triangles[index],
), 'order 6 render topology differs from native');
const order6XyzMetric = compareFloat(
  'order6.render_xyz', order6RenderXyz, order6Truth.fields.render_xyz,
);
const order6LogMetric = compareFloat(
  'order6.render_log_error', order6RenderLog, order6Truth.fields.render_log_error,
);
require(order6XyzMetric.maxAbs <= 1e-12,
        `order 6 render_xyz maxAbs=${order6XyzMetric.maxAbs}`);
// A few near-zero field values amplify native/WASM libm differences after
// log10.  Pair the pointwise ceiling with whole-field mean and RMS bounds.
require(order6LogMetric.maxAbs <= 2.0 &&
        order6LogMetric.meanAbs <= 2e-3 &&
        order6LogMetric.rms <= 3e-2 &&
        Math.abs(order6LogMetric.meanSigned - (-0.00039088793731142354)) <= 1e-8,
        `order 6 render_log_error metrics=${JSON.stringify(order6LogMetric)}`);
console.log('WASM_ORDER6_SMOKE_OK',
            `nsrc=${wasm._solver_result_nsrc()}`,
            `nrender=${wasm._solver_result_nrender()}`,
            `grf=${wasm._solver_result_grf_error()}`,
            `renderXyzMaxAbs=${order6XyzMetric.maxAbs}`,
            `renderLogMaxAbs=${order6LogMetric.maxAbs}`);
wasm._solver_clear();

callStage('_solver_run', 4n, 12n, 4n, 0n, 0.1,
          0n, 1e-3, ...refs4.pointers);
require(Number(wasm._solver_result_nsrc()) === 1360,
        `4x12 source count=${wasm._solver_result_nsrc()}, expected 1360`);
wasm._solver_clear();

require(solverRun(refs4, 12n, 36n, 4n, 2n, 0.1) === 107,
        'solver must reject an unknown surface id');
wasm._solver_clear();
require(solverRun(refs4, 12n, 36n, 4n, 1n, 0.0) === 108,
        'solver must reject a non-positive W7-X restol');
wasm._solver_clear();
require(solverRun(refs4, 12n, 36n, 4n, 1n, Number.NaN) === 108,
        'solver must reject a non-finite W7-X restol');
wasm._solver_clear();

let w7xDirectWasm;
if (process.env.WASM_SKIP_W7X !== '1') {
  callStage('_solver_run', 12n, 36n, 6n, 1n, 0.1,
            0n, 1e-3, ...refs6.pointers);
  const w7xNsrc = Number(wasm._solver_result_nsrc());
  const w7xNrender = Number(wasm._solver_result_nrender());
  const w7xNfaces = Number(wasm._solver_result_ntriangles());
  require(w7xNsrc === w7xTruth.nsrc,
          `w7x source count=${w7xNsrc}, expected ${w7xTruth.nsrc}`);
  require(w7xNrender === w7xTruth.nrender,
          `w7x render count=${w7xNrender}, expected ${w7xTruth.nrender}`);
  require(w7xNfaces === w7xTruth.nfaces,
          `w7x face count=${w7xNfaces}, expected ${w7xTruth.nfaces}`);
  require(Math.abs(wasm._solver_result_grf_error() - w7xTruth.grfError) < 1e-8,
          `w7x GRF=${wasm._solver_result_grf_error()} differs from native`);
  const w7xRenderXyz = copyF64('_solver_copy_render_xyz', 3 * w7xNrender);
  const w7xRenderLog = copyF64('_solver_copy_render_log_error', w7xNrender);
  const w7xTriangles = copyI64('_solver_copy_render_triangles', 3 * w7xNfaces);
  w7xDirectWasm = {
    grfError: wasm._solver_result_grf_error(),
    u: copyF64('_solver_copy_u', w7xNsrc),
    renderTriangles: w7xTriangles,
  };
  require(w7xTriangles.every(
    (value, index) => value === w7xTruth.fields.render_triangles[index],
  ), 'w7x render topology differs from native');
  const w7xXyzMetric = compareFloat(
    'w7x.render_xyz', w7xRenderXyz, w7xTruth.fields.render_xyz,
  );
  const w7xLogMetric = compareFloat(
    'w7x.render_log_error', w7xRenderLog, w7xTruth.fields.render_log_error,
  );
  require(w7xXyzMetric.maxAbs <= 1e-12,
          `w7x render_xyz maxAbs=${w7xXyzMetric.maxAbs}`);
  // One of 82,500 render vertices sits near the 1e-16 floor, where libm
  // noise swings the pointwise log by ~2 decades (measured 2.096).  The mean,
  // RMS, and pinned signed-mean signature carry the real sensitivity.
  require(w7xLogMetric.maxAbs <= 3.0 &&
          w7xLogMetric.meanAbs <= 2e-3 &&
          w7xLogMetric.rms <= 3e-2 &&
          Math.abs(w7xLogMetric.meanSigned - (-0.0000360063885458944)) <= 1e-8,
          `w7x render_log_error metrics=${JSON.stringify(w7xLogMetric)}`);
  console.log('WASM_W7X_SMOKE_OK',
              `nsrc=${w7xNsrc}`,
              `nrender=${w7xNrender}`,
              `grf=${wasm._solver_result_grf_error()}`,
              `renderXyzMaxAbs=${w7xXyzMetric.maxAbs}`,
              `renderLogMeanSigned=${w7xLogMetric.meanSigned}`);
  wasm._solver_clear();
}

function copyF64(name, count) {
  const bytes = count * 8;
  const pointer = wasm._malloc(bytes);
  require(pointer !== 0, `${name}: malloc failed`);
  try {
    const status = wasm[name](pointerArgument(pointer), BigInt(count));
    require(status === 0, `${name} failed: ${status}`);
    return new Float64Array(wasm.HEAPU8.buffer, pointer, count).slice();
  } finally {
    wasm._free(pointerArgument(pointer));
  }
}

function copyI64(name, count) {
  const bytes = count * 8;
  const pointer = wasm._malloc(bytes);
  require(pointer !== 0, `${name}: malloc failed`);
  try {
    const status = wasm[name](pointerArgument(pointer), BigInt(count));
    require(status === 0, `${name} failed: ${status}`);
    return new BigInt64Array(wasm.HEAPU8.buffer, pointer, count).slice();
  } finally {
    wasm._free(pointerArgument(pointer));
  }
}

function collect() {
  return collectFixture(truth);
}

function collectFixture(expected) {
  const nsrc = Number(wasm._solver_result_nsrc());
  const nrender = Number(wasm._solver_result_nrender());
  const nfaces = Number(wasm._solver_result_ntriangles());
  require(nsrc === expected.nsrc, `nsrc=${nsrc}, expected ${expected.nsrc}`);
  require(nrender === expected.nrender,
          `nrender=${nrender}, expected ${expected.nrender}`);
  require(nfaces === expected.nfaces, `nfaces=${nfaces}, expected ${expected.nfaces}`);
  return {
    grfError: wasm._solver_result_grf_error(),
    fields: {
      sx: copyF64('_solver_copy_sx', 3 * nsrc),
      snx: copyF64('_solver_copy_snx', 3 * nsrc),
      sw: copyF64('_solver_copy_sw', nsrc),
      ub: copyF64('_solver_copy_ub', nsrc),
      ubn: copyF64('_solver_copy_ubn', nsrc),
      u: copyF64('_solver_copy_u', nsrc),
      render_xyz: copyF64('_solver_copy_render_xyz', 3 * nrender),
      render_log_error: copyF64('_solver_copy_render_log_error', nrender),
      render_triangles: copyI64('_solver_copy_render_triangles', 3 * nfaces),
    },
  };
}

function runOnce() {
  wasm._solver_clear();
  callStage('_solver_run', 12n, 36n, 4n, 0n, 0.1,
            0n, 1e-3, ...refs4.pointers);
  return collect();
}

function compareFloat(name, actual, expected) {
  require(actual.length === expected.length, `${name}: length mismatch`);
  let maxAbs = 0, maxScaled = 0, referenceMax = 0;
  let sumAbs = 0, sumSquare = 0, sumSigned = 0;
  for (let i = 0; i < actual.length; ++i) {
    require(Number.isFinite(actual[i]), `${name}[${i}] is non-finite`);
    const abs = Math.abs(actual[i] - expected[i]);
    maxAbs = Math.max(maxAbs, abs);
    maxScaled = Math.max(maxScaled, abs / Math.max(1, Math.abs(expected[i])));
    referenceMax = Math.max(referenceMax, Math.abs(expected[i]));
    sumAbs += abs;
    sumSquare += abs * abs;
    sumSigned += actual[i] - expected[i];
  }
  return {
    maxAbs,
    maxScaled,
    fieldScaled: maxAbs / Math.max(1, referenceMax),
    meanAbs: sumAbs / actual.length,
    meanSigned: sumSigned / actual.length,
    rms: Math.sqrt(sumSquare / actual.length),
  };
}

function hashF64(values) {
  const bytes = Buffer.allocUnsafe(8 * values.length);
  for (let i = 0; i < values.length; ++i) bytes.writeDoubleLE(values[i], 8 * i);
  return createHash('sha256').update(bytes).digest('hex');
}

const tolerances = {
  sx: { maxAbs: 1e-12 },
  snx: { maxAbs: 1e-12 },
  sw: { maxAbs: 1e-12 },
  ub: { maxAbs: 1e-11 },
  ubn: { maxAbs: 1e-11 },
  u: { fieldScaled: 3e-3 },
  render_xyz: { maxAbs: 1e-12 },
  render_log_error: { maxAbs: 1.2 },
  grfError: { maxScaled: 1e-6 },
};

// The log field contains a few near-zero native values whose pointwise log can
// move by roughly one decade while the underlying field remains close. Its
// broad max bound is therefore paired with a deterministic signed-mean
// signature from the pinned compiler/toolchain. This makes the gate sensitive
// to a single corrupted render value instead of allowing a +1.0 perturbation.
const metricSignatures = {
  render_log_error: {
    meanSigned: { expected: 0.00001402499314912925, tolerance: 1e-9 },
  },
};

const deterministicFieldHashes = {
  u: 'ce4ae30ba941d9e62646cf93c838b202c1d6715753c3e90751824b0f10f83728',
  render_log_error: 'ae4e6605a7c7b348383fdf6300a3a915bd9235af714b7c87dc69909b818fc3f0',
};

function enforceTolerance(label, name, metric) {
  for (const [measure, limit] of Object.entries(tolerances[name])) {
    require(metric[measure] <= limit,
      `${label}: ${name}.${measure}=${metric[measure]} exceeds ${limit}`);
  }
  for (const [measure, signature] of Object.entries(metricSignatures[name] || {})) {
    require(Math.abs(metric[measure] - signature.expected) <= signature.tolerance,
      `${label}: ${name}.${measure}=${metric[measure]} differs from ` +
      `${signature.expected} by more than ${signature.tolerance}`);
  }
}

function compareRun(label, run) {
  const metrics = {};
  for (const name of Object.keys(truth.fields)) {
    if (name === 'render_triangles') {
      require(run.fields[name].every((value, index) => value === truth.fields[name][index]),
              `${label}: triangle topology differs`);
    } else {
      metrics[name] = compareFloat(name, run.fields[name], truth.fields[name]);
    }
  }
  metrics.grfError = compareFloat('grfError',
      new Float64Array([run.grfError]), new Float64Array([truth.grfError]));
  for (const [name, metric] of Object.entries(metrics)) {
    enforceTolerance(label, name, metric);
  }
  for (const [name, expected] of Object.entries(deterministicFieldHashes)) {
    require(hashF64(run.fields[name]) === expected,
      `${label}: ${name} deterministic SHA-256 differs`);
  }
  return metrics;
}

const first = runOnce();
const firstMemory = wasm.HEAPU8.buffer.byteLength;
if (process.env.WASM_PERTURB_RENDER_LOG === '1') {
  first.fields.render_log_error[0] += 1.0;
}
if (process.env.WASM_PERTURB_RENDER_LOG_BALANCED === '1') {
  first.fields.render_log_error[0] += 0.5;
  first.fields.render_log_error[1] -= 0.5;
}
if (process.env.WASM_PRINT_FIELD_HASHES === '1') {
  for (const name of Object.keys(deterministicFieldHashes)) {
    console.log('WASM_FIELD_SHA256', name, hashF64(first.fields[name]));
  }
  wasm._solver_clear();
  refs4.dispose();
  refs6.dispose();
  process.exit(0);
}
if (process.env.WASM_DIAGNOSTIC === '1') {
  for (const [name, values] of Object.entries(first.fields)) {
    if (name === 'render_triangles') continue;
    let nonFinite = 0, firstNonFinite = -1, min = Infinity, max = -Infinity;
    for (let i = 0; i < values.length; ++i) {
      if (!Number.isFinite(values[i])) {
        ++nonFinite;
        if (firstNonFinite < 0) firstNonFinite = i;
      } else {
        min = Math.min(min, values[i]);
        max = Math.max(max, values[i]);
      }
    }
    console.log('WASM_FIELD_DIAG', JSON.stringify({ name, nonFinite, firstNonFinite, min, max }));
  }
  console.log('WASM_GRF_DIAG', first.grfError);
  wasm._solver_clear();
  refs4.dispose();
  refs6.dispose();
  process.exit(0);
}
const firstMetrics = compareRun('run1', first);
if (process.env.WASM_SINGLE_RUN === '1') {
  console.log('WASM_PARITY_METRICS', JSON.stringify(firstMetrics));
  console.log('WASM_SINGLE_RUN_OK', `nsrc=${truth.nsrc}`, `nrender=${truth.nrender}`,
              `grf=${first.grfError}`);
  wasm._solver_clear();
  refs4.dispose();
  refs6.dispose();
  process.exit(0);
}

const second = runOnce();
const secondMemory = wasm.HEAPU8.buffer.byteLength;
const secondMetrics = compareRun('run2', second);
for (const name of Object.keys(first.fields)) {
  const a = first.fields[name], b = second.fields[name];
  require(a.length === b.length && a.every((value, index) => value === b[index]),
          `${name}: repeated WASM run differs`);
}
wasm._solver_clear();
callStage('_solver_run', 12n, 36n, 4n, 0n, 0.1,
          0n, 1e-3, ...refs4.pointers);
const convenience = collect();
const convenienceMemory = wasm.HEAPU8.buffer.byteLength;
compareRun('convenience', convenience);
if (process.env.WASM_SKIP_MEMORY_STABILITY !== '1') {
  require(convenienceMemory === secondMemory,
          `WASM memory still grows after warmup: ${firstMemory}, ${secondMemory}, ${convenienceMemory}`);
}

console.log('WASM_PARITY_METRICS', JSON.stringify(firstMetrics));
console.log('WASM_ORDER4_PARITY_OK', `nsrc=${truth.nsrc}`, `nrender=${truth.nrender}`,
            `grf=${first.grfError}`, `memory=${convenienceMemory}`);
wasm._solver_clear();

if (fmmFixtureDir) {
  const order4DirectWasm = {
    grfError: first.grfError,
    u: first.fields.u,
    renderTriangles: first.fields.render_triangles,
  };
  const fmmCases = [
    {
      label: 'builtin-order4-e3', file: 'builtin-order4-e3.bin',
      refs: refs4, args: [12n, 36n, 4n, 0n, 0.1, 1n, 1e-3],
      direct: truth, wasmDirect: order4DirectWasm, eps: 1e-3,
    },
    {
      label: 'builtin-order6-e6', file: 'builtin-order6-e6.bin',
      refs: refs6, args: [12n, 36n, 6n, 0n, 0.1, 1n, 1e-6],
      direct: order6Truth, wasmDirect: order6DirectWasm, eps: 1e-6,
    },
    {
      label: 'w7x-order6-e6', file: 'w7x-order6-e6.bin',
      refs: refs6, args: [12n, 36n, 6n, 1n, 0.1, 1n, 1e-6],
      direct: w7xTruth, wasmDirect: w7xDirectWasm, eps: 1e-6,
    },
  ];
  const selectedFmmCases = requestedFmmCases
    ? fmmCases.filter((testCase) => requestedFmmCases.has(testCase.label))
    : fmmCases;
  if (requestedFmmCases) {
    require(selectedFmmCases.length === requestedFmmCases.size,
      `unknown WASM_FMM_CASES value: ${process.env.WASM_FMM_CASES}`);
  }
  for (const testCase of selectedFmmCases) {
    const nativeFmm = parseFixture(path.join(fmmFixtureDir, testCase.file));
    const limit = 20 * testCase.eps;
    require(testCase.wasmDirect,
      `${testCase.label}: missing WASM Direct control result`);
    require(nativeFmm.ntri === testCase.direct.ntri &&
            nativeFmm.nsrc === testCase.direct.nsrc &&
            nativeFmm.nrender === testCase.direct.nrender &&
            nativeFmm.nfaces === testCase.direct.nfaces,
            `${testCase.label}: native FMM topology dimensions differ from Direct`);
    require(nativeFmm.fields.render_triangles.every(
      (value, index) => value === testCase.direct.fields.render_triangles[index]),
    `${testCase.label}: native FMM topology differs from Direct`);

    const nativeDirect = compareFloat(
      `${testCase.label}.native-direct.u`, nativeFmm.fields.u,
      testCase.direct.fields.u,
    );
    require(nativeDirect.fieldScaled <= limit,
      `${testCase.label}: native FMM/Direct u error ${nativeDirect.fieldScaled} > ${limit}`);
    require(Math.abs(nativeFmm.grfError-testCase.direct.grfError) <= limit,
      `${testCase.label}: native FMM/Direct GRF error exceeds ${limit}`);

    wasm._solver_clear();
    callStage('_solver_run', ...testCase.args, ...testCase.refs.pointers);
    const wasmFmm = collectFixture(nativeFmm);
    require(wasmFmm.fields.render_triangles.every(
      (value, index) => value === nativeFmm.fields.render_triangles[index]),
    `${testCase.label}: WASM FMM topology differs from native FMM`);
    require(wasmFmm.fields.render_triangles.every(
      (value, index) => value === testCase.wasmDirect.renderTriangles[index]),
    `${testCase.label}: WASM FMM topology differs from WASM Direct`);
    for (const [name, values] of Object.entries(wasmFmm.fields)) {
      if (name === 'render_triangles') continue;
      require(values.every(Number.isFinite),
        `${testCase.label}: WASM FMM ${name} contains non-finite values`);
    }
    const wasmNativeRaw = compareFloat(
      `${testCase.label}.wasm-native.u`, wasmFmm.fields.u, nativeFmm.fields.u,
    );
    const wasmDirect = compareFloat(
      `${testCase.label}.wasm-direct.u`, wasmFmm.fields.u,
      testCase.wasmDirect.u,
    );
    require(wasmDirect.fieldScaled <= limit,
      `${testCase.label}: WASM FMM/Direct u error ${wasmDirect.fieldScaled} > ${limit}`);
    require(Math.abs(wasmFmm.grfError-testCase.wasmDirect.grfError) <= limit,
      `${testCase.label}: WASM FMM/Direct GRF error exceeds ${limit}`);

    const nativeEffect = new Float64Array(nativeFmm.nsrc);
    const wasmEffect = new Float64Array(nativeFmm.nsrc);
    for (let i = 0; i < nativeFmm.nsrc; ++i) {
      nativeEffect[i] = nativeFmm.fields.u[i]-testCase.direct.fields.u[i];
      wasmEffect[i] = wasmFmm.fields.u[i]-testCase.wasmDirect.u[i];
    }
    const effectMetric = compareFloat(
      `${testCase.label}.fmm-effect`, wasmEffect, nativeEffect,
    );
    require(effectMetric.fieldScaled <= limit,
      `${testCase.label}: cross-toolchain FMM effect error ` +
      `${effectMetric.fieldScaled} > ${limit}`);
    console.log('WASM_FMM_PARITY_OK', testCase.label,
                `ntri=${nativeFmm.ntri}`,
                `nsrc=${nativeFmm.nsrc}`,
                `uNativeDirect=${nativeDirect.fieldScaled}`,
                `uWasmDirect=${wasmDirect.fieldScaled}`,
                `effect=${effectMetric.fieldScaled}`,
                `rawCrossToolchain=${wasmNativeRaw.fieldScaled}`,
                `nativeGrf=${nativeFmm.grfError}`,
                `wasmGrf=${wasmFmm.grfError}`,
                `heap=${wasm.HEAPU8.buffer.byteLength}`);
    wasm._solver_clear();
  }
}

// --- Progress-event characterization (Task 7) -------------------------------
// Capture the raw progress events emitted during a real built-in order-4 solve
// and prove their semantics: monotone loop counters, exactly one 100% per long
// phase, the exact stage ordering, and a single result event whose value is the
// solver's own returned GRF.  This is a durable characterization of the final
// WASM artifact, not a re-derivation of any numerical quantity.
function captureProgress(refs, runArgs) {
  const events = [];
  wasm.onSolverProgress = (stage, current, total, aux0, aux1, value) => {
    events.push({
      stage: Number(stage),
      current: Number(current),
      total: Number(total),
      aux0: Number(aux0),
      aux1: Number(aux1),
      value,
    });
  };
  try {
    const status = solverRun(refs, ...runArgs);
    require(status === 0,
      `_solver_run failed: status=${status} last=${wasm._solver_last_error()}`);
  } finally {
    delete wasm.onSolverProgress;
  }
  return events;
}

function stageEvents(events, stage) {
  return events.filter((event) => event.stage === stage);
}

function requireOnce(events, stage) {
  const matches = stageEvents(events, stage);
  require(matches.length === 1, `expected exactly one stage ${stage}, got ${matches.length}`);
  return matches[0];
}

// Long loops (direct blocks, close-count panels, RRQ panels): monotone counters,
// bounds 0 < current <= total, and exactly one 100% event.
function requireLongPhase(events, progressStage, label) {
  const progress = stageEvents(events, progressStage);
  require(progress.length >= 1, `${label}: no progress events`);
  let previous = 0;
  let final = 0;
  for (const event of progress) {
    require(event.current > 0 && event.current <= event.total,
            `${label}: counter ${event.current}/${event.total} out of bounds`);
    require(event.current > previous,
            `${label}: counter not strictly increasing (${previous} -> ${event.current})`);
    previous = event.current;
    if (event.current === event.total) ++final;
  }
  require(final === 1, `${label}: expected exactly one 100% event, got ${final}`);
  return progress[progress.length - 1];
}

const progressEvents = captureProgress(refs4, [12n, 36n, 4n, 0n, 0.1]);
const progressGrf = wasm._solver_result_grf_error();
const progressNsrc = Number(wasm._solver_result_nsrc());

const geometryBegin = requireOnce(progressEvents, 1);
const geometryReady = requireOnce(progressEvents, 2);
require(geometryReady.total === truth.nsrc,
        `geometry-ready nsrc=${geometryReady.total}, expected ${truth.nsrc}`);
const directBegin = requireOnce(progressEvents, 3);
const directFinal = requireLongPhase(progressEvents, 4, 'direct');
require(directBegin.aux0 === 0 && directFinal.aux0 === 0,
        'Direct progress events must carry kernel id 0');
const countBegin = requireOnce(progressEvents, 5);
const countFinal = requireLongPhase(progressEvents, 6, 'close/count');
const rrqBegin = requireOnce(progressEvents, 7);
const rrqFinal = requireLongPhase(progressEvents, 8, 'close/rrq');
const scatter = requireOnce(progressEvents, 9);
const render = requireOnce(progressEvents, 10);
const result = requireOnce(progressEvents, 11);

// Exact stage ordering by position in the emission stream.
const orderedMilestones = [
  ['geometry begin', geometryBegin],
  ['geometry ready', geometryReady],
  ['direct begin', directBegin],
  ['direct 100%', directFinal],
  ['close count begin', countBegin],
  ['close count 100%', countFinal],
  ['close RRQ begin', rrqBegin],
  ['close RRQ 100%', rrqFinal],
  ['scatter', scatter],
  ['render', render],
  ['result', result],
];
for (let i = 1; i < orderedMilestones.length; ++i) {
  const [prevName, prevEvent] = orderedMilestones[i - 1];
  const [name, event] = orderedMilestones[i];
  const prevIndex = progressEvents.indexOf(prevEvent);
  const index = progressEvents.indexOf(event);
  require(prevIndex < index, `stage order violated: ${prevName} must precede ${name}`);
}

require(result.value === progressGrf,
        `result event value=${result.value} differs from _solver_result_grf_error()=${progressGrf}`);
require(result.current === progressNsrc,
        `result event current=${result.current}, expected nsrc=${progressNsrc}`);
wasm._solver_clear();

const fmmProgressEvents = captureProgress(
  refs4, [12n, 36n, 4n, 0n, 0.1, 1n, 1e-3]);
const fmmBegin = requireOnce(fmmProgressEvents, 3);
const fmmComplete = requireOnce(fmmProgressEvents, 4);
require(fmmBegin.aux0 === 1 && fmmComplete.aux0 === 1,
        'FMM progress events must carry kernel id 1');
require(fmmComplete.current === 1 && fmmComplete.total === 1,
        `FMM completion must be 1/1, got ${fmmComplete.current}/${fmmComplete.total}`);
require(Number.isFinite(fmmComplete.value) && fmmComplete.value >= 0,
        `FMM completion time is invalid: ${fmmComplete.value}`);
const fmmFollowingStages = [5, 6, 7, 8, 9, 10, 11];
let previousFmmIndex = fmmProgressEvents.indexOf(fmmComplete);
for (const stage of fmmFollowingStages) {
  const event = stage === 6 || stage === 8
    ? requireLongPhase(fmmProgressEvents, stage, `FMM stage ${stage}`)
    : requireOnce(fmmProgressEvents, stage);
  const index = fmmProgressEvents.indexOf(event);
  require(previousFmmIndex < index,
          `FMM completion must precede stage ${stage}`);
  previousFmmIndex = index;
}
console.log('WASM_FMM_PROGRESS_OK',
            `events=${fmmProgressEvents.length}`,
            `seconds=${fmmComplete.value}`);
wasm._solver_clear();

callStage('_solver_run', 12n, 36n, 4n, 0n, 0.1,
          1n, 1e-3, ...refs4.pointers);
const fmmMemoryAfterSecond = wasm.HEAPU8.buffer.byteLength;
wasm._solver_clear();
callStage('_solver_run', 12n, 36n, 4n, 0n, 0.1,
          1n, 1e-3, ...refs4.pointers);
const fmmMemoryAfterThird = wasm.HEAPU8.buffer.byteLength;
require(fmmMemoryAfterThird === fmmMemoryAfterSecond,
  `WASM FMM memory still grows after warmup: ` +
  `${fmmMemoryAfterSecond}, ${fmmMemoryAfterThird}`);
console.log('WASM_FMM_MEMORY_OK', `memory=${fmmMemoryAfterThird}`);
wasm._solver_clear();

// A callback that throws on every event must not alter solver status or results.
wasm.onSolverProgress = () => { throw new Error('progress handler intentionally throws'); };
let throwingStatus;
try {
  throwingStatus = solverRun(refs4, 12n, 36n, 4n, 0n, 0.1);
} finally {
  delete wasm.onSolverProgress;
}
require(throwingStatus === 0,
        `throwing progress handler changed solver status to ${throwingStatus}`);
const throwingRun = collect();
compareRun('throwing-progress', throwingRun);
wasm._solver_clear();

console.log('WASM_PROGRESS_ORDER_OK',
            `events=${progressEvents.length}`,
            `directBlocks=${directFinal.total}`,
            `ntri=${countFinal.total}`,
            `grf=${result.value}`);
refs4.dispose();
refs6.dispose();
