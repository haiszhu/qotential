function [A11 A12 A13 A22 A23 A33] = Sto3dSLPnmat_mex(t,s)

r0 = t.x;
r0n = t.nx;
m = numel(r0(1,:));
r = s.x;
n = numel(r(1,:));
w = s.w;

mex_id_ = 'sto3dslpnmat(c i int[x], c i double[xx], c i double[xx], c i int[x], c i double[xx], c i double[x], c o double[xx], c o double[xx], c o double[xx], c o double[xx], c o double[xx], c o double[xx])';
[A11, A12, A13, A22, A23, A33] = specialquad(mex_id_, m, r0, r0n, n, r, w, 1, 3, m, 3, m, 1, 3, n, n, m, n, m, n, m, n, m, n, m, n, m, n);

end

