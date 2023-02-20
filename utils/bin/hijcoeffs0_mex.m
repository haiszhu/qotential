function [coeffs_ext,pow_ext,len] = hijcoeffs0_mex(i,j,order2)

mex_id_ = 'hijcoeffs0(i int[x], i int[x], i int[x], o double[x], o int[xx], o int[x])';
[coeffs_ext, pow_ext, len] = specialquad(mex_id_, i, j, order2, 1, 1, 1, order2, 3, order2, 1);
coeffs_ext = coeffs_ext';

end

