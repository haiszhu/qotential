function coeffs_updt = hijcoeffsk_mex(i, j, k, pow)

kp1 = k+1;
mex_id_ = 'lenofnckcoeffs(i int[x], i int[x], i int[x], o int[x], o int[x])';
[len, idx_flag] = specialquad(mex_id_, i, j, k, 1, 1, 1, 1, kp1);

[yPow,xPow]=meshgrid(0:j,0:i);
pow = [xPow(:),yPow(:)]'; 

ip1 = i+1; jp1 = j+1;
ip1xjp1 = ip1*jp1; 

mex_id_ = 'hijcoeffsk0(i int[x], i int[x], i int[x], i int[xx], i int[x], o double[x])';
[coeffs_updt] = specialquad(mex_id_, i, j, k, pow, len, 1, 1, 1, 2, ip1xjp1, 1, ip1xjp1);
coeffs_updt = coeffs_updt';
end

