function [Amu, mat_perm, Q_mat_perm, Amu_c] = Sto3dModifiedMu(r,nr,convention_flag)
% modified Stokes SLP eval matrix (to make kernel exact form)
% only once per patch
% 
% [A, mat_perm] = Sto3dModifiedmu(r,nr)
%
% modified density mu_tilde = ( A(:,1:end/4)...
%                             + A(:,end/4+1:end/2)*r0(1)...
%                             + A(:,end/2+1:3*end/4)*r0(2)...
%                             + A(:,3*end/4+1:end)*r0(3)...
%                                ) * mat_perm * [mu_i mu_j mu_k]'
% out of this has ordering: [4 entry; 4; ...]
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
%  Q_mat_perm - (4*N * 4*N) density map from ordering [4;4;4;...;4] to [N;N;N;N] 
% 
% Hai 08/24/22 (only modify 1st row for now, i.e. 1st component in velocity)
% 
% Hai 08/26/22 allow just output Amu that follows vectored kernel convention

if nargin < 3, convention_flag = 0; end

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

%   Amu0(k,k:n:end) = Mat_10{k}(1,:);
%   Amu0(n+k,k:n:end) = Mat_10{k}(2,:);
%   Amu0(2*n+k,k:n:end) = Mat_10{k}(3,:);
%   Amu0(3*n+k,k:n:end) = Mat_10{k}(4,:);
end
Amu = [blkdiag(Mat_10{:}), blkdiag(Mat_11{:}), blkdiag(Mat_12{:}), blkdiag(Mat_13{:})]; % 
Amu = -Amu/2; % modified density (extra -1/2 from stokes slp mu to laplace dlp mu_tilde,...
              % due to kernel coefficient 1/(4*pi) and 1/(8*pi) and quaternion -1 for scalar component after multiplication)

% create map from density ordering of [mu1(:);mu2(:);mu3(:)] to ordering of [mu1(1);mu2(1);mu3(1); ... mu1(n);mu2(n);mu3(n)]
idx1 = 1:3*n+3:(n-1)*(3*n+3)+1; % nonzeros elements index to map mu1(:) to spaced by 3 version...
idx2 = idx1 + n*3*n + 1;        % nonzeros elements index to map mu2(:) to spaced by 3 version...
idx3 = idx2 + n*3*n + 1;        % nonzeros elements index to map mu3(:) to spaced by 3 version...
row1 = mod(idx1-1,3*n)+1; row2 = mod(idx2-1,3*n)+1; row3 = mod(idx3-1,3*n)+1;
mat_perm = sparse([row1 row2 row3],1:3*n,ones(1,3*n)); 
% mat_perm = sparse(3*n,3*n); mat_perm(idx1) = 1; mat_perm(idx2) = 1; mat_perm(idx3) = 1; % equivalently, above is faster
% Mu = mat_perm*[mu_i mu_j mu_k]'; % check...

% create map from mu_tilde ordering of [mu_tilde1(1);mu_tilde2(1);mu_tilde3(1);mu_tilde4(1); ... mu_tilde1(n);mu_tilde2(n);mu_tilde3(n);mu_tilde4(n)]
% to ordering of [mu_tilde1(:);mu_tilde2(:);mu_tilde3(:);mu_tilde4(:);
idx0 = 1:4*n+4:(n-1)*(4*n+4)+1; 
idx1 = idx0 + n*4*n + 1;
idx2 = idx1 + n*4*n + 1;
idx3 = idx2 + n*4*n + 1;
row0 = mod(idx0-1,4*n)+1; row1 = mod(idx1-1,4*n)+1; row2 = mod(idx2-1,4*n)+1; row3 = mod(idx3-1,4*n)+1;
Q_mat_perm = sparse([row0 row1 row2 row3],1:4*n,ones(1,4*n))'; 

% instead, output Amu that respects vectored kernel convention
if convention_flag
  Amu0 = Q_mat_perm*Amu(:,1:end/4)*mat_perm; Amu1 = Q_mat_perm*Amu(:,end/4+1:end/2)*mat_perm; 
  Amu2 = Q_mat_perm*Amu(:,end/2+1:3*end/4)*mat_perm; Amu3 = Q_mat_perm*Amu(:,3*end/4+1:end)*mat_perm; 
  Amu_c = [Amu0 Amu1 Amu2 Amu3];
else
  Amu_c = [];
end

end