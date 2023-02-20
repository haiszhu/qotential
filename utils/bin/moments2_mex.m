function [Nk Mk Lk]=moments2_mex(r0,r,order)
% single-threaded MEX

% MATLAB preprocessing here
N = numel(r)/3;
flag = double(nargout); % crucial even though it's int type in Fortran
                        % later for Lk output?
% since Nk Mk has type 'output' of given size, no MATLAB preallocation is needed
% Code for calling the Fortran function
orderp1 = order + 1; r0 = r0'; r = r';
mex_id_ = 'moments2(i double[x], i int, i double[xx], i int, i int, o double[xx], o double[xx], o double[xx])';
[Nk, Mk, Lk] = specialquad(mex_id_, r0, N, r, order, flag, 3, N, 3, N, orderp1, N, orderp1, N, orderp1);

end


% @function [Omega_x,Omega_y,Omega_z,Omega21,Omega22,Omega23,Omega31,Omega32,Omega33] = omegaallDLPn_mex(r0,M_all,L_all,order,onm_0_rs_x,onm_1_rs_x,onm_2_rs_x,onm_3_rs_x,onm_0_rs_y,onm_1_rs_y,onm_2_rs_y,onm_3_rs_y,onm_0_rs_z,onm_1_rs_z,onm_2_rs_z,onm_3_rs_z,h_dim,ijIdx)
% 
% m = numel(r0(1,:));
% M_all = reshape(M_all,[],m);
% L_all = reshape(L_all,[],m);
% dim1 = numel(M_all(:,1));
% n = numel(onm_0_rs_x(:,1))/h_dim;
% nh_dim = n*h_dim;
% morder = 2*order + 2;
% h_dim4 = h_dim*4;
% 
% if nargout <=3
% # FORTRAN omegaalldlpn(int[1] m, double[3,m] r0, int[1] dim1, double[dim1,m] M_all, double[dim1,m] L_all, int[1] n, int[1] morder, double[nh_dim,11] onm_0_rs_x, double[nh_dim,11] onm_1_rs_x, double[nh_dim,11] onm_2_rs_x, double[nh_dim,11] onm_3_rs_x, double[nh_dim,11] onm_0_rs_y, double[nh_dim,11] onm_1_rs_y, double[nh_dim,11] onm_2_rs_y, double[nh_dim,11] onm_3_rs_y, double[nh_dim,11] onm_0_rs_z, double[nh_dim,11] onm_1_rs_z, double[nh_dim,11] onm_2_rs_z, double[nh_dim,11] onm_3_rs_z, int[1] h_dim, int[2,h_dim] ijIdx, output double[m,h_dim4] Omega_x, output double[m,h_dim4] Omega_y, output double[m,h_dim4] Omega_z);
% else
% # FORTRAN omegaalldlpnmat(int[1] m, double[3,m] r0, int[1] dim1, double[dim1,m] M_all, double[dim1,m] L_all, int[1] n, int[1] morder, double[nh_dim,11] onm_0_rs_x, double[nh_dim,11] onm_1_rs_x, double[nh_dim,11] onm_2_rs_x, double[nh_dim,11] onm_3_rs_x, double[nh_dim,11] onm_0_rs_y, double[nh_dim,11] onm_1_rs_y, double[nh_dim,11] onm_2_rs_y, double[nh_dim,11] onm_3_rs_y, double[nh_dim,11] onm_0_rs_z, double[nh_dim,11] onm_1_rs_z, double[nh_dim,11] onm_2_rs_z, double[nh_dim,11] onm_3_rs_z, int[1] h_dim, int[2,h_dim] ijIdx, output double[m,h_dim4] Omega_x, output double[m,h_dim4] Omega_y, output double[m,h_dim4] Omega_z, output double[m,h_dim4] Omega21, output double[m,h_dim4] Omega22, output double[m,h_dim4] Omega23, output double[m,h_dim4] Omega31, output double[m,h_dim4] Omega32, output double[m,h_dim4] Omega33);
% end
% 
% end

