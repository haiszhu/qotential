function [coeffs,pow] = hijcoeffs_mex(i,j)

ip1 = i+1; jp1 = j+1;
ip1xjp1 = ip1*jp1;
mex_id_ = 'hijcoeffs(i int[x], i int[x], o double[x], o int[xx])';
[coeffs, pow] = specialquad(mex_id_, i, j, 1, 1, ip1xjp1, 3, ip1xjp1);
coeffs = coeffs';

end

