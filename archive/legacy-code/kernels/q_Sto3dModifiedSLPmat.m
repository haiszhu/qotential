function [q_A1, q_mu] = q_Sto3dModifiedSLPmat(t,s,mu)
% dense modified quaternion Stokes SLP Nystrom eval matrix from s.x to t.x
%
% q_A = q_Sto3dModifiedSLPmat(t,s,sigma)
%
% q_A1 = 1/(8*pi*|r0-r|^3)*(r0-r)*nr*conj(nr)*...
%           ( conj(r0-r)*sigma(1) - sigma*(r0(1)-r(1)) ) 
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, nx (3*N) normals, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes
%  mu - density (3*N)
% Outputs:
%  q_A - (3*M*N) matrix getting potentials from density values
%       right now only 1st entry, so matrix size M*N 
%  q_mu - (1*N) quaternionic density
%
% Hai 08/25/22

% quaternion
q_r0 = quaternion(zeros(1,size(t.x,2)),t.x(1,:),t.x(2,:),t.x(3,:));
q_r  = quaternion(zeros(1,size(s.x,2)),s.x(1,:),s.x(2,:),s.x(3,:));
q_nr = quaternion(zeros(1,size(s.nx,2)),s.nx(1,:),s.nx(2,:),s.nx(3,:));
q_mu = quaternion(zeros(1,size(mu,2)),mu(1,:),mu(2,:),mu(3,:));

% form matrix
q_A1 = 1/(8*pi)*(q_r0(:).*q_nr-q_r.*q_nr)./norm(q_r0(:)-q_r).^3.*...
      ( ((conj(q_nr(:)).*conj(q_r0)).'-conj(q_nr).*conj(q_r)).*(s.w.*mu(1,:)./norm(q_mu).^2.*conj(q_mu))...
        - (((conj(q_nr(:)).*s.w(:)).*t.x(1,:)).'-conj(q_nr).*s.x(1,:).*s.w) );


end