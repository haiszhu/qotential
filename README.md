# qotential

tools for high‑order quaternionic product quadratures that resolve singular and nearly singular 3D Laplace layer potentials. The implementation follows Ref. 1 and Ref. 2.

## Layout
- `kernels/`: dense Laplace/Stokes kernel evaluation (inspired and partialy copied from https://github.com/ahbarnett/BIE3D), and FMM wrapper with similar interface as the dense kernel evaluation, quaternionic evaluation variants.
- `utils/`: a few smooth quadrature utilities, harmonic bases, FMM (https://github.com/flatironinstitute/FMM3D), kd‑tree (https://github.com/taiya/kdtree).
- `utils/f`: source files. `lap3ddlp.f` implements the original close‑evaluation scheme in Ref. 1; `lap3dslp.f` implements the correct Laplace SLP formulation in Ref. 2. The improved RRQ scheme is in `src/`.
- `src/`: Fortran implementation of the RRQ close evaluation of Ref. 2 — per‑patch uniform refinement, the S^{lm} tc integrals, and the assembled `rrq_r64`, with a MATLAB gateway in `matlab/`.
- `external/`: [LineQuaaadrature](https://github.com/haiszhu/LineQuaaadrature) (line quadrature, ssq weights, solid angle) and [QuatApproximation](https://github.com/haiszhu/QuatApproximation) (S^{lm} moments, translation operators), as submodules. Build both before `make lib`.
- `test/`: `stellarator_bench.f90` (throughput), `stellarator_grf.m` (built-in stellarator GRF accuracy), and `w7x_grf.m` (W7-X). Their surface meshes are generated in memory by the LineQuaaadrature Fortran/MEX geometry interfaces.

## References
1. Hai Zhu, and Shravan Veerapaneni. 2022. “High-Order Close Evaluation of Laplace Layer Potentials: A Differential Geometric Approach.” *SIAM Journal on Scientific Computing*.
2. Shidong Jiang, and Hai Zhu. 2024. “Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions.” *arXiv preprint arXiv:2411.08342*.

## To do

* computation pipeline and interface redesign
* solid angle fix
* quaternion approximation on curved element (thm)
* Helmholtz kernel
* Stokes kernel
* quaternion operation abstraction?
