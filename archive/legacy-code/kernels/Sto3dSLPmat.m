function A = Sto3dSLPmat(t,s)
% dense Stokes SLP Nystrom eval matrix from sources to targets
%
% A = Sto3dSLPmat(t,s)
%
% Inputs:
%  s - source surface struct with fields:
%      x (3*N) nodes, w (1*N) quadr weights
%  t - target struct with fields:
%      x (3*M) nodes
% Outputs:
%  A - (3*M*3*N) matrix getting potentials from density values
%
% NB: fills 3*M*3*N matrices, so assumes small cases. 
%

% Hai 08/23/22

x=t.x; y=s.x; M = numel(x)/3;
d1 = bsxfun(@minus,x(1,:)',y(1,:));   % 3 coords of displacement matrix (M*N)
d2 = bsxfun(@minus,x(2,:)',y(2,:));
d3 = bsxfun(@minus,x(3,:)',y(3,:));
ir = 1./sqrt(d1.^2+d2.^2+d3.^2);   % 1/r matrix, size M*N

cross12 = d1.*d2; cross13 = d1.*d3; cross23 = d2.*d3; 
A = [ir + d1.^2.*ir.^3,    cross12.*ir.^3,    cross13.*ir.^3;...
        cross12.*ir.^3, ir + d2.^2.*ir.^3,    cross23.*ir.^3;...
        cross13.*ir.^3,    cross23.*ir.^3, ir + d3.^2.*ir.^3];
A = 1/(8*pi)*A .* repmat([s.w(:)' s.w(:)' s.w(:)'], [3*M 1]);    

end