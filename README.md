# rrq

tools for high‑order quaternionic product quadratures that resolve singular and nearly singular 3D Laplace layer potentials. The implementation follows Ref. 1 and Ref. 2.

## Layout
- `src/` and `matlab/`: maintained Fortran RRQ implementation and MATLAB MEX gateway.
- `utils/`: maintained MATLAB quadrature helpers.
- `external/`: [LineQuaaadrature](https://github.com/haiszhu/LineQuaaadrature), [QuatApproximation](https://github.com/haiszhu/QuatApproximation), and [kd-tree](https://github.com/taiya/kdtree) submodules.
- `test/`: maintained native and MATLAB stellarator tests, with surface meshes generated in memory by the LineQuaaadrature Fortran/MEX geometry interfaces.
- `web/stellarator/`: browser solver and its build and test pipeline.
- `archive/legacy-code/`: inactive historical MATLAB demos and kernel helpers; this directory is not added by `setup.m`.

## Build and benchmark

```bash
git submodule update --init --recursive
make lib
make -C test
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ./build/stellarator_bench 5 128 0
```

This runs orders 4–12 with the 128-wide SIMD kernel and uniform refinement.

The torus comparison additionally requires a built FMM3DBIE checkout.

```bash
make fmm3dbie-mex FMM3DBIE_DIR=~/git/fmm3dbie
```

Compile the platform-specific kdtree MEX files once from MATLAB:

```matlab
cd('/path/to/qotential/external/kdtree/toolbox')
kdtree_compile
```

## References
1. Hai Zhu, and Shravan Veerapaneni. 2022. “High-Order Close Evaluation of Laplace Layer Potentials: A Differential Geometric Approach.” *SIAM Journal on Scientific Computing*.
2. Shidong Jiang, and Hai Zhu. 2024. “Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions.” *arXiv preprint arXiv:2411.08342*.

## To do list

* (to do) check examples on different systems, and update readme
* (to do) minimal Helmholtz slp + dlp
* (to do) minimal Stokes slp...
* (to do) check, and move all line quadrature to LQ submodule, move all quaternion approximation and locloc related to QA
* (to do) solid angle fix
* (to do) quaternion approximation on curved element (thm)
* (to do) laplace slpn, dlpn, maybe (no legacy code available...) slpnn, dlpnn 
* (to do) helmholtz slpn + dlpn
* (to do) stokes (no legacy code available...) slpn, dlp, dlpn...