function [nck_coeffs,nck_pow] = polylapder2d_mex(i, j, ell)

nck_der = 2*[ell:-1:0; 0:1:ell];
idx = (nck_der(1,:)<=i) & (nck_der(2,:)<=j);
len = sum(idx);
nck_der = nck_der(:,idx);
nck_idx = 1:ell+1; nck_idx = nck_idx(idx);
nck_pow = [];
ellp1 = ell+1;
mex_id_ = 'lenofnckcoeffs(i int[x], i int[x], i int[x], o int[x], o int[x])';
[len, idx_flag] = specialquad(mex_id_, i, j, ell, 1, 1, 1, 1, ellp1);

if len
mex_id_ = 'polylapder2d(i int[x], i int[x], i int[x], i int[x], i int[xx], i int[x], o double[x], o int[xx])';
[nck_coeffs, nck_pow] = specialquad(mex_id_, i, j, ell, len, nck_der, nck_idx, 1, 1, 1, 1, 2, len, len, len, 3, len);
else
nck_coeffs = []; nck_pow = [];
end
end

