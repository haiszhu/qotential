function A = Lap3dDLPmat_mex(t,s)

r0 = t.x;
m = numel(r0(1,:));
r = s.x;
rn = s.nx;
n = numel(r(1,:));
w = s.w;

mex_id_ = 'lap3ddlpmat(i int[x], i double[xx], i int[x], i double[xx], i double[xx], i double[x], o double[xx])';
[A] = specialquad(mex_id_, m, r0, n, r, rn, w, 1, 3, m, 1, 3, n, 3, n, n, m, n);

end


