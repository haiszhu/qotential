function [qnm_i_1,qnm_i_2,qnm_i_3] = qnm_i(r,gradF,i,lptype,gradxyz)
% q^{(nm)}_i: (qnm_i_1,qnm_i_2,qnm_i_3) vector
% gradF = \nabla \nu
% wrapper for previously separately qnm function, follow notation from
% Shidong's notes, without kernel denominator part (moments come later)
% 
% same for SLPn (lptype 't') & DLP (lptype 'd'), 
% 
% f1*\nu*g1 + f2*\nu*g2 for SLP (lptype 's')
%
% Hai 11/27/22, add support for SLP, Hai 11/29/22
% add support for DLPn, Hai 01/15/23, need gradxyz (partial x, y, or z)

if nargin < 4, lptype = 'd'; end

if lptype == 'd' || lptype == 't'
  switch i  % DLP & SLPn
      case 0
          [qnm_i_1,qnm_i_2,qnm_i_3] = q0kl(struct('x',r),gradF,'d');
      case 1
          [qnm_i_1,qnm_i_2,qnm_i_3] = q1kl(struct('x',r),gradF,'d');
      case 2
          [qnm_i_1,qnm_i_2,qnm_i_3] = q2kl(struct('x',r),gradF,'d');
      case 3
          [qnm_i_1,qnm_i_2,qnm_i_3] = q3kl(struct('x',r),gradF,'d');
  end
elseif lptype == 's' % SLP
  switch i
      case 0
          [qnm_i_1,qnm_i_2,qnm_i_3] = q0kl(struct('x',r),gradF,'s');
      case 1
          [qnm_i_1,qnm_i_2,qnm_i_3] = q1kl(struct('x',r),gradF,'s');
      case 2
          [qnm_i_1,qnm_i_2,qnm_i_3] = q2kl(struct('x',r),gradF,'s');
      case 3
          [qnm_i_1,qnm_i_2,qnm_i_3] = q3kl(struct('x',r),gradF,'s');
  end
elseif lptype == 'T' % DLPn
  switch i
    case 0
      [qnm_i_1,qnm_i_2,qnm_i_3] = q0kl_DLPn(struct('x',r),gradF,gradxyz);
    case 1
      [qnm_i_1,qnm_i_2,qnm_i_3] = q1kl_DLPn(struct('x',r),gradF,gradxyz);
    case 2
      [qnm_i_1,qnm_i_2,qnm_i_3] = q2kl_DLPn(struct('x',r),gradF,gradxyz);
    case 3
      [qnm_i_1,qnm_i_2,qnm_i_3] = q3kl_DLPn(struct('x',r),gradF,gradxyz);
  end
end

end