import path from 'node:path';
import { pathToFileURL } from 'node:url';

const modulePath = process.argv[2];
if (!modulePath) throw new Error('usage: memory64-pointer.runner.mjs <module.mjs>');
const imported = await import(pathToFileURL(path.resolve(modulePath)).href);
const wasm = await imported.default();

const offset = wasm._malloc(16);
if (typeof offset !== 'number' || !Number.isSafeInteger(offset) || offset <= 0) {
  throw new Error(`_malloc did not return a numeric heap offset: ${String(offset)}`);
}
const rawPointer = BigInt(offset);
const returned = wasm._pointer_identity(rawPointer);
if (typeof returned !== 'bigint' || returned !== rawPointer) {
  throw new Error(`raw pointer ABI mismatch: ${typeof returned} ${String(returned)}`);
}
wasm._pointer_store(rawPointer, 0x5a);
if (wasm.HEAPU8[offset] !== 0x5a) {
  throw new Error('wasm64 pointer did not address the expected heap byte');
}

function checkedOffset(pointer) {
  if (typeof pointer !== 'bigint' || pointer < 0n ||
      pointer > BigInt(Number.MAX_SAFE_INTEGER)) {
    throw new RangeError('wasm64 pointer is outside the safe JavaScript offset range');
  }
  return Number(pointer);
}
if (checkedOffset(returned) !== offset) throw new Error('checked pointer conversion failed');
try {
  checkedOffset(BigInt(Number.MAX_SAFE_INTEGER) + 1n);
  throw new Error('unsafe pointer conversion was accepted');
} catch (error) {
  if (!(error instanceof RangeError)) throw error;
}

wasm._free(rawPointer);
console.log('MEMORY64_POINTER_ABI_OK',
  `malloc=${typeof offset}`, `raw=${typeof returned}`, `offset=${offset}`);
