# qotential

tools for high‑order quaternionic product quadratures that resolve singular and nearly singular 3D layer potentials (Laplace, .. ). The implementation follows Ref. 1 and Ref. 2.

## Layout
- `kernels/`: dense Laplace/Stokes kernel evaluation (inspired and partialy copied from https://github.com/ahbarnett/BIE3D), and FMM wrapper with similar interface as the dense kernel evaluation, quaternionic evaluation variants.
- `utils/`: a few smooth quadrature utilities, harmonic bases, FMM (https://github.com/flatironinstitute/FMM3D), kd‑tree (https://github.com/taiya/kdtree).
- `utils/f`: source files. `lap3ddlp.f` implements the original close‑evaluation scheme in Ref. 1; `lap3dslp.f` implements the correct Laplace SLP formulation in Ref. 2, but not the improved RRQ scheme...

## References
1. Hai Zhu, and Shravan Veerapaneni. 2022. “High-Order Close Evaluation of Laplace Layer Potentials: A Differential Geometric Approach.” *SIAM Journal on Scientific Computing*.
2. Shidong Jiang, and Hai Zhu. 2024. “Recursive reduction quadrature for the evaluation of Laplace layer potentials in three dimensions.” *arXiv preprint arXiv:2411.08342*.

Please consider citing the above works if this code benefits your research.
