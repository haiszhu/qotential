function [nck_coeffs,nck_pow] = polylapder2d_mex(i, j, ell)

nck_der = 2*[ell:-1:0; 0:1:ell];
idx = (nck_der(1,:)<=i) & (nck_der(2,:)<=j);
len = sum(idx);
nck_der = nck_der(:,idx);
nck_idx = 1:ell+1; nck_idx = nck_idx(idx);
nck_pow = [];
ellp1 = ell+1;
mex_id_ = 'lenofnckcoeffs(c i int[x], c i int[x], c i int[x], c o int[x], c o int[x])';
[len, idx_flag] = specialquad(mex_id_, i, j, ell, 1, 1, 1, 1, ellp1);

if len
mex_id_ = 'polylapder2d(c i int[x], c i int[x], c i int[x], c i int[x], c i int[xx], c i int[x], c o double[x], c o int[xx])';
[nck_coeffs, nck_pow] = specialquad(mex_id_, i, j, ell, len, nck_der, nck_idx, 1, 1, 1, 1, 2, len, len, len, 3, len);
else
nck_coeffs = []; nck_pow = [];
end
end

