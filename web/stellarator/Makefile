LF ?= lfortran
EMCC ?= emcc
FMM3D_SRC ?= $(HOME)/git/FMM3D

.PHONY: wasm fmm3d-wasm fmm3d-test test parity app sanitize serve

fmm3d-wasm:
	EMCC="$(EMCC)" scripts/build-fmm3d-wasm.sh

fmm3d-test:
	bash test/fmm3d-manifest.test.sh
	bash test/fort2c-allocate-stat.test.sh
	bash test/fmm3d-build-cache.test.sh

wasm:
	LF="$(LF)" EMCC="$(EMCC)" scripts/build-wasm.sh

test:
	npm test

parity: wasm
	$(MAKE) -C ../../test sgrf_fmm FMM3D_DIR="$(FMM3D_SRC)"
	@tmp_dir=$$(mktemp -d "$${TMPDIR:-/tmp}/stellarator-fmm-parity.XXXXXX"); \
	  trap 'rm -rf "$$tmp_dir"' EXIT; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../build/stellarator_fmm \
	    "$$tmp_dir/builtin-order4-e3.bin" 12 36 4 builtin 0.1 1e-3; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../build/stellarator_fmm \
	    "$$tmp_dir/builtin-order6-e6.bin" 12 36 6 builtin 0.1 1e-6; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../build/stellarator_fmm \
	    "$$tmp_dir/w7x-order6-e6.bin" 12 36 6 w7x 0.1 1e-6; \
	  python3 test/native-fmm-parity.py --fixture-dir "$$tmp_dir"; \
	  WASM_FMM_FIXTURE_DIR="$$tmp_dir" node test/wasm-api-node.mjs

app: wasm test
	npm run build

sanitize:
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O1 WASM_SANITIZE=address scripts/build-wasm.sh
	WASM_SINGLE_RUN=1 node test/wasm-api-node.mjs
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O2 scripts/build-wasm.sh

serve: wasm
	npm run dev
