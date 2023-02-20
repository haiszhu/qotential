function [q0kl_i,q0kl_j,q0kl_k] = q0kl(s,gradF,lptype)
% quaternion representation of q_0^(k,l)
%     vectored form: f^(k,l)(r) cross product (r0-r)
%     q_0^(k,l)_x part (vectored function)
%              _y part
%              _z part
%     _x part, _y part, _z part parts go into quaternion
%     each part, say _x part, consists of fixed, x0-dependent, y0-dependent
%     , z0-dependent parts, therefore the last dimension is 4
%     to make it target independent: dim(q0) = N * num of coeffs * 4
%     to compute q_0^(k,l) = q0(:,:,1) + q0(:,:,2)*tx(1) + q0(:,:,3)*tx(2) + q0(:,:,4)*tx(3);
%                q_0^(k,l) is quaternion 
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes
%  gradF - struct containing gradient of harmonic basis
%          gradF.F1 (N * num of coeffs) partial derivative of basis w.r.t. x 
%          gradF.F2 ... w.r.t. y
%          gradF.F3 ... w.r.t. z
%
% Hai 08/30/22
% no longer quaternion

if nargin < 3, lptype = 'd'; end % by default, DLP/SLPn

switch lptype
  case 'd'
    F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    tmp = zeros(size(F1)); tmp2 = zeros(size(F1,1),size(F1,2),4);
    q0kl_i = tmp2; q0kl_j = tmp2; q0kl_k = tmp2;
    q0kl_i(:,:,1) = (-s.x(2,:)').*F3 - (-s.x(3,:)').*F2;
    q0kl_i(:,:,2) = tmp; q0kl_i(:,:,3) = F3; q0kl_i(:,:,4) = -F2;
    q0kl_j(:,:,1) = (-s.x(3,:)').*F1 - (-s.x(1,:)').*F3;
    q0kl_j(:,:,2) = - F3; q0kl_j(:,:,3) = tmp; q0kl_j(:,:,4) = F1;
    q0kl_k(:,:,1) = (-s.x(1,:)').*F2 - (-s.x(2,:)').*F1;
    q0kl_k(:,:,2) = F2; q0kl_k(:,:,3) = - F1; q0kl_k(:,:,4) = tmp;
  case 's' % require harmonic polynomial in addition to its gradient
    F = gradF.F; F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
    tmp = zeros(size(F1)); tmp2 = zeros(size(F1,1),size(F1,2),5); 
    q0kl_i = tmp2; q0kl_j = tmp2; q0kl_k = tmp2;
    q0kl_i(:,:,1) = s.x(1,:)'.*F; q0kl_i(:,:,2) = -F;  q0kl_i(:,:,3) = tmp; q0kl_i(:,:,4) = tmp;
    q0kl_i(:,:,5) = F1;
    q0kl_j(:,:,1) = s.x(2,:)'.*F; q0kl_j(:,:,2) = tmp; q0kl_j(:,:,3) = -F;  q0kl_j(:,:,4) = tmp;
    q0kl_j(:,:,5) = F2;
    q0kl_k(:,:,1) = s.x(3,:)'.*F; q0kl_k(:,:,2) = tmp; q0kl_k(:,:,3) = tmp; q0kl_k(:,:,4) = -F;
    q0kl_k(:,:,5) = F3;
end

end