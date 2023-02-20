function [Omega_x,Omega_y,Omega_z] = omegaallvec_mex(r0,M_all,order,onm_0_rs,onm_1_rs,onm_2_rs,onm_3_rs,h_dim,ijIdx)

m = numel(r0(1,:));
M_all = reshape(M_all,[],m);
dim1 = numel(M_all(:,1));
n = numel(onm_0_rs(:,1))/h_dim;
nh_dim = n*h_dim;
morder = 2*order + 2;
h_dim4 = h_dim*4;

mex_id_ = 'omegaallvec(i int[x], i double[xx], i int[x], i double[xx], i int[x], i int[x], i double[xx], i double[xx], i double[xx], i double[xx], i int[x], i int[xx], o double[xx], o double[xx], o double[xx])';
[Omega_x, Omega_y, Omega_z] = specialquad(mex_id_, m, r0, dim1, M_all, n, morder, onm_0_rs, onm_1_rs, onm_2_rs, onm_3_rs, h_dim, ijIdx, 1, 3, m, 1, dim1, m, 1, 1, nh_dim, 4, nh_dim, 4, nh_dim, 4, nh_dim, 4, 1, 2, h_dim, m, h_dim4, m, h_dim4, m, h_dim4);

end

