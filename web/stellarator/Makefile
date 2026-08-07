LF ?= lfortran
EMCC ?= emcc

.PHONY: wasm test parity app sanitize serve

wasm:
	LF="$(LF)" EMCC="$(EMCC)" scripts/build-wasm.sh

test:
	npm test

parity: wasm
	node test/wasm-api-node.mjs

app: wasm test
	npm run build

sanitize:
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O1 WASM_SANITIZE=address scripts/build-wasm.sh
	WASM_SINGLE_RUN=1 node test/wasm-api-node.mjs
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O2 scripts/build-wasm.sh

serve: wasm
	npm run dev
