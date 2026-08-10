function [q_A, q_A_basis, A_basis] = q_Lap3dDLPmat(t,s,gradF)
% dense modified quaternion Laplace DLP Nystrom eval matrix from s.x to t.x
%
% q_A = q_Lap3dDLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, nx (3*N) unit normals, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes
%  gradF - struct containing gradient of harmonic basis
%          gradF.F1 (N*basis order) partial derivative w.r.t. x 
%
% Outputs:
%  q_A - (M*N) matrix getting potentials from density values
%
%  q_A_basis - A times basis, for j-th target q_A_basis(:,:,j)
%
%  A_basis - projected on dx\wedge dy, dy\wedge dz, dz\wedge dx ... (somewhat more complicated)
%           q0, three components of f(r) \cross (r0-r) in quaternion form
%           q1, three components of -((r0-r)\dot f^(k,l))e_1 + (x0-x)*f^(k,l) + f_1^(k,l)(r0-r) in quaternion form
%           q2, three components of -((r0-r)\dot f^(k,l))e_2 + (y0-y)*f^(k,l) + f_2^(k,l)(r0-r) in quaternion form
% 
% Hai 08/25/22, add quaternion kernel times basis for debugging, Hai
% 08/30/22

if nargin==0, test_q_Lap3dDLPmat; return; end

% quaternion
q_r0 = quaternion(zeros(1,size(t.x,2)),t.x(1,:),t.x(2,:),t.x(3,:));
q_r  = quaternion(zeros(1,size(s.x,2)),s.x(1,:),s.x(2,:),s.x(3,:));
q_nr = quaternion(zeros(1,size(s.nx,2)),s.nx(1,:),s.nx(2,:),s.nx(3,:));

% dist
r0mr_norm3 = norm(q_r0(:)-q_r).^3;

% form matrix
q_A = -1/(4*pi)*(q_r0(:).*q_nr.*s.w-q_r.*q_nr.*s.w)./r0mr_norm3;

% form matrix times basis (still matrix or tensor?)
% for verification purpose, close eval is target dependent, so maybe of
% dimension N * basis order * M
if nargin > 2
  F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3; F0 = zeros(size(F1));
  q_basis = quaternion(F0,F1,F2,F3);
  tmp = zeros(size(F1,1),size(F1,2),size(t.x,2));
  q_A_basis = quaternion(tmp,tmp,tmp,tmp);
  for j=1:size(t.x,2) % j-th target
    q_A_basis(:,:,j) = q_A(j,:).'.*q_basis;
  end

  A_basis = [];
  Q_i = zeros(size(F1,1),size(F1,2),size(t.x,2));
  Q_j = zeros(size(F1,1),size(F1,2),size(t.x,2));
  Q_k = zeros(size(F1,1),size(F1,2),size(t.x,2));
  % A_basis_0: three components of alpha_0^(k,l)
  for j=1:size(t.x,2) % j-th target
    % f^(k,l)(r) cross product (r0-r)
    txj = t.x(:,j);
    Q_i(:,:,j) = (txj(2)-s.x(2,:)').*F3 - (txj(3)-s.x(3,:)').*F2;
    Q_j(:,:,j) = (txj(3)-s.x(3,:)').*F1 - (txj(1)-s.x(1,:)').*F3;
    Q_k(:,:,j) = (txj(1)-s.x(1,:)').*F2 - (txj(2)-s.x(2,:)').*F1;
  end
  A_basis.q0 = quaternion(tmp,Q_i,Q_j,Q_k);
  % A_basis_1: three components of alpha_1^(k,l)
  for j=1:size(t.x,2) % j-th target
    % -((r0-r)\dot f^(k,l))e_1 + (x0-x)*f^(k,l) + f_1^(k,l)(r0-r)
    txj = t.x(:,j);
    r0mr_dot_f = (txj(1)-s.x(1,:)').*F1 + (txj(2)-s.x(2,:)').*F2 + (txj(3)-s.x(3,:)').*F3;
    Q_i(:,:,j) = -(-r0mr_dot_f + 2*(txj(1)-s.x(1,:)').*F1);
    Q_j(:,:,j) = -((txj(1)-s.x(1,:)').*F2 + (txj(2)-s.x(2,:)').*F1);
    Q_k(:,:,j) = -((txj(1)-s.x(1,:)').*F3 + (txj(3)-s.x(3,:)').*F1);
  end
  A_basis.q1 = quaternion(tmp,Q_i,Q_j,Q_k);
  % A_basis_2: three components of alpha_2^(k,l)
  for j=1:size(t.x,2) % j-th target
    % -((r0-r)\dot f^(k,l))e_2 + (y0-y)*f^(k,l) + f_2^(k,l)(r0-r)
    txj = t.x(:,j);
    r0mr_dot_f = (txj(1)-s.x(1,:)').*F1 + (txj(2)-s.x(2,:)').*F2 + (txj(3)-s.x(3,:)').*F3;
    Q_i(:,:,j) = -((txj(2)-s.x(2,:)').*F1 + (txj(1)-s.x(1,:)').*F2);
    Q_j(:,:,j) = -(-r0mr_dot_f + 2*(txj(2)-s.x(2,:)').*F2);
    Q_k(:,:,j) = -((txj(2)-s.x(2,:)').*F3 + (txj(3)-s.x(3,:)').*F2);
  end
  A_basis.q2 = quaternion(tmp,Q_i,Q_j,Q_k);
  % A_basis_3: three components of alpha_3^(k,l)
  for j=1:size(t.x,2) % j-th target
    % -((r0-r)\dot f^(k,l))e_3 + (z0-z)*f^(k,l) + f_3^(k,l)(r0-r)
    txj = t.x(:,j);
    r0mr_dot_f = (txj(1)-s.x(1,:)').*F1 + (txj(2)-s.x(2,:)').*F2 + (txj(3)-s.x(3,:)').*F3;
    Q_i(:,:,j) = -((txj(3)-s.x(3,:)').*F1 + (txj(1)-s.x(1,:)').*F3);
    Q_j(:,:,j) = -((txj(3)-s.x(3,:)').*F2 + (txj(2)-s.x(2,:)').*F3);
    Q_k(:,:,j) = -(-r0mr_dot_f + 2*(txj(3)-s.x(3,:)').*F3);
  end
  A_basis.q3 = quaternion(tmp,Q_i,Q_j,Q_k);
  % keyboard
  A_basis.dist3 = r0mr_norm3;
  A_basis.w = s.w;

else
  q_A_basis = []; A_basis = [];
end

end

%%%%%%%%%
function test_q_Lap3dDLPmat
% math test: tau=-1 GRF for torus surface...
type = 'torus', so.a = 1; so.b = 0.5; o = [];
so.mp = 8; so.np = 12;
[pan N] = create_panels(type,so,o);
[r nr w] = getallnodes(pan);
s.x = r; s.nx = nr; s.w = w; s.N = N;

xin = [0.9; -0.2; 0.1]; xout = [1.9; 0.7; 1.0];    % both "far" from surf
t.x = [xin,xout]; t.nx = randn(3,2);               % target pointset, rand n
tau = -ones(s.N,1);           % DLP density
q_tau = quaternion(tau',zeros(1,numel(tau)),zeros(1,numel(tau)),zeros(1,numel(tau)));
q_A = q_Lap3dDLPmat(t,s);
q_u = q_A.*q_tau;
Q_u1 = compact(q_u(1,:)); Q_u2 = compact(q_u(2,:));
fprintf('torus DLP[-1] errs (val,grad): int %.3g, ext %.3g\n',sum(Q_u1(:,1))-1,sum(Q_u2(:,1)))

disp('timing test...')
ns = 5e3;
nt = 1e4;
s = []; s.x = randn(3,ns); s.w = ones(1,ns);
s.nx = randn(3,ns); s.nx = bsxfun(@times,s.nx,1./sqrt(sum(s.nx.^2,1)));
%norm(s.nx(:,1))-1
t.x = randn(3,nt); t.nx = randn(3,nt);
tic;
%profile clear; profile on
A = Lap3dDLPmat(t,s);
%[A An] = Lap3dDLPmat(t,s);   % twice as slow, understandably
%profile off; profile viewer
t1=toc;
q_A = q_Lap3dDLPmat(t,s);
t2 = toc-t1;
fprintf('filled lap dipole pot mat with quaternion is %d as slow\n',t2/t1)
end