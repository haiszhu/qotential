function [q1kl_i,q1kl_j,q1kl_k] = q1kl(s,gradF,lptype)
% quaternion representation of q_1^(k,l)
%     vectored form: 
%     to make it target independent: dim(q0) = N * num of coeffs * 4
%     to compute q_1^(k,l) = q1(:,:,1) + q1(:,:,2)*tx(1) + q1(:,:,3)*tx(2) + q1(:,:,4)*tx(3);
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes
%  gradF - struct containing gradient of harmonic basis
%          gradF.F1 (N * num of coeffs) partial derivative of basis w.r.t. x 
%          gradF.F2 ... w.r.t. y
%          gradF.F3 ... w.r.t. z
%
% Hai 08/30/33

if nargin < 3, lptype = 'd'; end % by default, DLP/SLPn

switch lptype
  case 'd'
    F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    tmp = zeros(size(F1)); tmp2 = zeros(size(F1,1),size(F1,2),4);
    q1kl_i = tmp2; q1kl_j = tmp2; q1kl_k = tmp2;
    mr_dot_f = (-s.x(1,:)').*F1 + (-s.x(2,:)').*F2 + (-s.x(3,:)').*F3;
    q1kl_i(:,:,1) = -(-mr_dot_f + 2*(-s.x(1,:)').*F1);
    q1kl_i(:,:,2) = -F1; q1kl_i(:,:,3) = F2; q1kl_i(:,:,4) = F3;
    q1kl_j(:,:,1) = -((-s.x(1,:)').*F2 + (-s.x(2,:)').*F1);
    q1kl_j(:,:,2) = -F2; q1kl_j(:,:,3) = -F1; q1kl_j(:,:,4) = tmp;
    q1kl_k(:,:,1) = -((-s.x(1,:)').*F3 + (-s.x(3,:)').*F1);
    q1kl_k(:,:,2) = -F3; q1kl_k(:,:,3) = tmp; q1kl_k(:,:,4) = -F1;
  case 's'
    F = gradF.F; F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    tmp = zeros(size(F1)); tmp2 = zeros(size(F1,1),size(F1,2),5); 
    q1kl_i = tmp2; q1kl_j = tmp2; q1kl_k = tmp2;
    q1kl_i; % 0
    q1kl_j(:,:,1) =  s.x(3,:)'.*F; q1kl_j(:,:,2) = tmp; q1kl_j(:,:,3) = tmp; q1kl_j(:,:,4) = -F;
    q1kl_j(:,:,5) = -F3;
    q1kl_k(:,:,1) = -s.x(2,:)'.*F; q1kl_k(:,:,2) = tmp; q1kl_k(:,:,3) = F;   q1kl_k(:,:,4) = tmp;
    q1kl_k(:,:,5) = F2;
end

end