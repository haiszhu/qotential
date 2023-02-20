function t_val = interpval(t,s_val,x0)
% using 1d interpolation, assume GL nodes.
% interpolate 3d coordinates s_val = (x,y,z) to subpanel parametrized by t
%
% easier to debug if interpolation is accurate for tiny patches
% assume RCIP works with 1d interpolation, then this should also work?
%
% interpolate from 3xN source to 3*Nt targets
% s_val is source value, assume x otimes x grid, where x is 1d GL
% t_val is target value, assume [t(1,:);t(2,:)]...
%
% 12/22/22 Hai

% (x,y,z) of s_val
p = length(x0);

% intptmp1 & intptmp2
intptmp = interpmat_1d([t(1,:)';t(2,:)'],x0); % 1d GL nodes along t(2) direction

% 1st interpolate to (t(1),x0), where x0 is gauss legendre nodes
% intptmp1 = interpmat_1d(t(1,:)',x0); % 1d GL nodes along t(2) direction
intptmp1 = intptmp(1:end/2,:);

% then interpolate to (t(1),t(2))
% intptmp2 = interpmat_1d(t(2,:)',x0); 
intptmp2 = intptmp(end/2+1:end,:);

% more general
t_val = zeros(size(s_val,1),size(t,2));
for k=1:size(s_val,1) % assume each row of s_val needs to be interpilated to target t
  s_valk = reshape(s_val(k,:),p,p);
  s_valtmp = intptmp1*s_valk';
  t_val(k,:) = sum((s_valtmp.*intptmp2)');
end

end

function L = interpmat_1d(t,s)
% INTERPMAT_1D   interpolation matrix from nodes in 1D to target nodes
%
% L = interpmat_1d(t,s) returns interpolation matrix taking values on nodes s
%  to target nodes t. Computed in Helsing style.

% bits taken from qplp/interpmatrix.m from 2013
% Barnett 7/17/16

if nargin==0, test_interpmat_1d; return; end

p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
end 