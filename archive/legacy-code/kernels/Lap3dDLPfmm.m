function [u,un,unx,uny,unz] = Lap3dDLPfmm(t,s,tau,fmm_eps)
% LAP3DDLPFMM. Laplace DLP fmm from sources to targets 
% similar interface as Lap3dDLPmat & Lap3dDLP_closepanel

% 
srcinfo.sources = s.x;
srcinfo.dipoles = tau'.*s.w.*s.nx; % srcinfo = rmfield(srcinfo,'charges'); if necessary
targ = t.x;
if nargin < 4, fmm_eps = 1e-15; end
ifppreg = 0;

if nargout < 2
  ifppregtarg = 1; % potential at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg/(4*pi))'; 
elseif nargout < 3 % need t.nx...
  ifppregtarg = 2; % potential and its gradient at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg/(4*pi))'; 
  un = (sum(t.nx.*U.gradtarg)/(4*pi))'; 
else
  ifppregtarg = 2; % potential and its gradient at target
  U = lfmm3d(fmm_eps,srcinfo,ifppreg,targ,ifppregtarg);
  u = (U.pottarg/(4*pi))'; 
  un = (sum(t.nx.*U.gradtarg)/(4*pi))'; 
  unx = (U.gradtarg(1,:)/(4*pi))';
  uny = (U.gradtarg(2,:)/(4*pi))';
  unz = (U.gradtarg(3,:)/(4*pi))';
end

end