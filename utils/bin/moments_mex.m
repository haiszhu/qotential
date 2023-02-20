function [Nk Mk]=moments_mex(r0,r,order)
% single-threaded MEX

% MATLAB preprocessing here
N = numel(r)/3;
flag = double(nargout); % crucial even though it's int type in Fortran
                        % later for Lk output?
% since Nk Mk has type 'output' of given size, no MATLAB preallocation is needed
% Code for calling the Fortran function
orderp1 = order + 1; r0 = r0'; r = r';
mex_id_ = 'moments(i double[x], i int, i double[xx], i int, i int, o double[xx], o double[xx])';
[Nk, Mk] = specialquad(mex_id_, r0, N, r, order, flag, 3, N, 3, N, orderp1, N, orderp1);

end

