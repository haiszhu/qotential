function [q2kl_i,q2kl_j,q2kl_k] = q2kl(s,gradF,lptype)
% quaternion representation of q_2^(k,l)
%     vectored form: 
%     to make it target independent: dim(q0) = N * num of coeffs * 4
%     to compute q_2^(k,l) = q1(:,:,1) + q1(:,:,2)*tx(1) + q1(:,:,3)*tx(2) + q1(:,:,4)*tx(3);
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
    q2kl_i = tmp2; q2kl_j = tmp2; q2kl_k = tmp2;
    mr_dot_f = (-s.x(1,:)').*F1 + (-s.x(2,:)').*F2 + (-s.x(3,:)').*F3;
    q2kl_i(:,:,1) = -((-s.x(2,:)').*F1 + (-s.x(1,:)').*F2);
    q2kl_i(:,:,2) = -F2; q2kl_i(:,:,3) = -F1; q2kl_i(:,:,4) = tmp;
    q2kl_j(:,:,1) = -(-mr_dot_f + 2*(-s.x(2,:)').*F2);
    q2kl_j(:,:,2) = F1; q2kl_j(:,:,3) = -F2; q2kl_j(:,:,4) = F3;
    q2kl_k(:,:,1) = -((-s.x(2,:)').*F3 + (-s.x(3,:)').*F2);
    q2kl_k(:,:,2) = tmp; q2kl_k(:,:,3) = -F3; q2kl_k(:,:,4) = -F2;
  case 's'
    F = gradF.F; F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    tmp = zeros(size(F1)); tmp2 = zeros(size(F1,1),size(F1,2),5); 
    q2kl_i = tmp2; q2kl_j = tmp2; q2kl_k = tmp2;
    q2kl_i(:,:,1) = -s.x(3,:)'.*F; q2kl_i(:,:,2) = tmp; q2kl_i(:,:,3) = tmp; q2kl_i(:,:,4) = F;
    q2kl_i(:,:,5) = F3;
    q2kl_j; % 0
    q2kl_k(:,:,1) =  s.x(1,:)'.*F; q2kl_k(:,:,2) = -F;  q2kl_k(:,:,3) = tmp; q2kl_k(:,:,4) = tmp;
    q2kl_k(:,:,5) = -F1;
end

end