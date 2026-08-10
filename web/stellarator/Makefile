LF ?= lfortran
EMCC ?= emcc
FMM3D_SRC ?= $(HOME)/git/FMM3D

.PHONY: wasm wasm64 artifacts fmm3d-wasm fmm3d-test memory64-test test parity parity64 parity64-long app sanitize serve

fmm3d-wasm:
	EMCC="$(EMCC)" scripts/build-fmm3d-wasm.sh

fmm3d-test:
	bash test/fmm3d-manifest.test.sh
	bash test/fort2c-allocate-stat.test.sh
	bash test/fort2c-real-complex-kind.test.sh
	bash test/fmm3d-build-cache.test.sh

memory64-test:
	bash test/memory64-pointer.test.sh

wasm:
	LF="$(LF)" EMCC="$(EMCC)" WASM_MEMORY64=0 scripts/build-wasm.sh

wasm64:
	LF="$(LF)" EMCC="$(EMCC)" WASM_MEMORY64=1 scripts/build-wasm.sh
	bash test/lfortran-wasm64-float128-symbols.test.sh

artifacts: wasm wasm64
	bash test/wasm-artifacts.test.sh
	scripts/write-wasm-provenance.sh
	scripts/update-wasm-checksums.sh
	cd public/wasm && shasum -a 256 --check SHA256SUMS

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

parity64: wasm64
	$(MAKE) -C ../../test sgrf_fmm FMM3D_DIR="$(FMM3D_SRC)"
	@tmp_dir=$$(mktemp -d "$${TMPDIR:-/tmp}/stellarator-fmm-parity64.XXXXXX"); \
	  trap 'rm -rf "$$tmp_dir"' EXIT; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../build/stellarator_fmm \
	    "$$tmp_dir/builtin-order6-e6.bin" 12 36 6 builtin 0.1 1e-6; \
	  python3 test/native-fmm-parity.py --fixture-dir "$$tmp_dir" \
	    --case builtin-order6-e6.bin; \
	  WASM_ARTIFACT=solver64 WASM_POINTER_BITS=64 WASM_SKIP_W7X=1 \
	    WASM_FMM_FIXTURE_DIR="$$tmp_dir" WASM_FMM_CASES=builtin-order6-e6 \
	    node test/wasm-api-node.mjs

parity64-long: wasm64
	$(MAKE) -C ../../test sgrf_fmm FMM3D_DIR="$(FMM3D_SRC)"
	@tmp_dir=$$(mktemp -d "$${TMPDIR:-/tmp}/stellarator-fmm-parity64-long.XXXXXX"); \
	  trap 'rm -rf "$$tmp_dir"' EXIT; \
	  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ../../build/stellarator_fmm \
	    "$$tmp_dir/builtin-order12-e12.bin" 48 144 12 builtin 0.1 1e-12; \
	  node test/wasm64-order12-parity.mjs \
	    "$$tmp_dir/builtin-order12-e12.bin"

app: artifacts test
	npm run build

sanitize:
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O1 WASM_SANITIZE=address scripts/build-wasm.sh
	WASM_SINGLE_RUN=1 node test/wasm-api-node.mjs
	LF="$(LF)" EMCC="$(EMCC)" WASM_OPT=-O2 scripts/build-wasm.sh

serve: artifacts
	npm run dev
