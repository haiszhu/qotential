function [A, mat_perm] = Sto3dModifiedSigma(r,nr)
% modified Stokes SLP eval matrix (to make kernel exact form)
% only once per patch
% 
% [A, mat_perm] = Sto3dModifiedSigma(r,nr)
%
% Inputs: 
%  r - source nodes (3*N), likely N = p^2 or p*(p+1)/2
%  nr - source normal (3*N)
% Outputs:
%  A - (4*N * 3*N*4) matrix getting quaternionic density values from
%       vector density values: 
%       4*N rows for quaternion representation (ordering: [4 entry;4;...])
%       3*N*1 columns for 3*N original density (ordering: [3 entry;3;...])
%       3*N*1 columns times r0(1) for 3*N original density, r0 - target
%       3*N*1 columns times r0(2) 
%       3*N*1 columns times r0(3)  
%  mat_perm - (3*N * 3*N) density map from ordering [N;N;N] to [3;3;3;...;3] 
% 
% Hai 08/24/22 (only modify 1st row for now, i.e. 1st component in velocity)
%

n = numel(r(:))/3;
% check modification matrix source by source
Mat_10 = cell(1,n); Mat_11 = cell(1,n); Mat_12 = cell(1,n); Mat_13 = cell(1,n);
mat_11 = [0,0,0;-2,0,0;0,-1,0;0,0,-1];   % to be multiplied by r0(1), target
mat_12 = [0,0,0;0,0,0;-1,0,0;0,0,0];     % r0(2)
mat_13 = [0,0,0;0,0,0;0,0,0;-1,0,0];     % r0(3)
for k=1:n
  mat_10 = [0,0,0;2*r(1,k),0,0;r(2,k),r(1,k),0;r(3,k),0,r(1,k)];
  nrmat = [0,nr(1,k),nr(2,k),nr(3,k);-nr(1,k),0,nr(3,k),-nr(2,k);-nr(2,k),-nr(3,k),0,nr(1,k);-nr(3,k),nr(2,k),-nr(1,k),0];
  % 
  Mat_10{k} = sparse(nrmat*mat_10);
  Mat_11{k} = sparse(nrmat*mat_11);
  Mat_12{k} = sparse(nrmat*mat_12);
  Mat_13{k} = sparse(nrmat*mat_13);
end
A = [blkdiag(Mat_10{:}), blkdiag(Mat_11{:}), blkdiag(Mat_12{:}), blkdiag(Mat_13{:})]; % 

% create map from density ordering of [sigma1(:);sigma2(:);sigma3(:)] to ordering of [sigma1(1);sigma2(1);sigma3(1); ... sigma1(n);sigma2(n);sigma3(n)]
idx1 = 1:3*n+3:(n-1)*(3*n+3)+1; % nonzeros elements index to map sigma1(:) to spaced by 3 version...
idx2 = idx1 + n*3*n + 1;        % nonzeros elements index to map sigma2(:) to spaced by 3 version...
idx3 = idx2 + n*3*n + 1;        % nonzeros elements index to map sigma3(:) to spaced by 3 version...
row1 = mod(idx1-1,3*n)+1; row2 = mod(idx2-1,3*n)+1; row3 = mod(idx3-1,3*n)+1;
mat_perm = sparse([row1 row2 row3],1:3*n,ones(1,3*n)); 
% mat_perm = sparse(3*n,3*n); mat_perm(idx1) = 1; mat_perm(idx2) = 1; mat_perm(idx3) = 1; % equivalently, above is faster
% Mu = mat_perm*[mu_i mu_j mu_k]'; % check...

end