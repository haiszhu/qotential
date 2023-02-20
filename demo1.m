% plain no-edge-adaptive close eval demo
% from 1 benign patch to close targets
%
% Hai 02/20/23
%

setup()
profile clear
profile on

% order
p = 8; % try <= 14, before condition number kicks in
ref = 1; % 1, 2, 3, 4...

% a patch on a torus
a = 1; b = 0.5; 
phi = outfun(a,b);
chart = @(t,p) [b*sin(t); (a+b*cos(t)).*sin(p); (a+b*cos(t)).*cos(p)]; 
orig = [0;0]; tpansiz1 = 1/2/ref; tpansiz2 = 1/4/ref;
[x0, w0, D] = gauss(p); [x1, x2] = meshgrid(x0); 
r = chart(orig(1)+tpansiz1*x1(:)',orig(2)+tpansiz2*x2(:)'); 
rc = mean(r,2); r = r - rc;
s = get_high_order_quad(r,p);
sf = get_high_order_quad(r,50*p); % upsampled 

% define close targets (num 27419) & density function
tmp = sqrt(sum(s.w))*(rand(3,100000)-1/2);
t.x = mean(s.x,2) + [1.25/3*tmp(1:2,:);1.25*tmp(3,:)]; % random targets
t.x = t.x + rc; t.x = t.x(:,phi(t.x(1,:),t.x(2,:),t.x(3,:))) - rc;
fmu = @(x,y,z) sin(3*x)+cos(2*y)+sin(z+1/2)+x.^2+y.^3+z.^4+exp(x.*y.*z);

%%% naive quadrature error
% naive quadr
A = Lap3dDLPmat_mex(t,s);
u = A*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
% upsampled quadr (this is 10 digit reference soln)
uf = Lap3dDLPfmm(t,sf,fmu(sf.x(1,:),sf.x(2,:),sf.x(3,:))',1e-13);
% error
err = abs(u-uf)/max(abs(uf(:)));
figure(1),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("naive eval: max abs err " + max(err) + " ")

%%% close evaluation quadrature error
% close eval quadr
[Ac,r] = Lap3dDLP_closepanel_demo(t,s,0,ref);
% Ac = Lap3dDLP_closepanel_demo(t,s,1);
% Ac = Lap3dDLP_closepanel(t,s,'e');
uc = Ac*fmu(s.x(1,:),s.x(2,:),s.x(3,:))';  
% error
err2 = abs(uc-uf)/max(abs(uf(:)));
figure(1), plot3(r(1,:),r(2,:),r(3,:),'.k-')
figure(2),clf,scatter3(t.x(1,:),t.x(2,:),t.x(3,:),7,log10(err2),'filled')
axis equal, hold on
colorbar, caxis([-10 0])
plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r')
title("close eval: max abs err " + max(err2) + " ")
plot3(r(1,:),r(2,:),r(3,:),'.k-')

profile viewer

keyboard

function [A,r] = Lap3dDLP_closepanel_demo(t,s,if_adapt,ref)
% close eval matrix, everything on-the-fly
%
% 
% Hai 02/20/23 

if nargin < 3, if_adapt = 0; end
if nargin < 4, ref = 1; end

% assume tensor product grid
order = sqrt(numel(s.x(1,:)));

% geometric process: shift, rescale, set origin, etc...
% compute length of pank "roughtly", and overwrite pank.x
sbd = []; sbd.p = ref*order; sbd.tpan = linspace(0,2*pi,9);
sbd.x_uvs_qntype = 'T';
sbd.xvr = s.x; 
sbd = quadr_3dline(sbd, [], 'G');
sx = s.x; %
r0 = t.x; %
r = sbd.x;
rtau = sbd.tang;

% set up interpolation matrix
[gradF,~] = evaltensorproductharmonicgrad_mex(sx,order);
h_dim = order^2;
n = numel(sx(:))/3;
Amu_c = [speye(n);sparse(3*n,n)];
F0 = zeros(h_dim); F1 = gradF.F1; F2 = gradF.F2; F3 = gradF.F3;
Mmatrix = [[ F0 -F1 -F2 -F3];...
           [ F1  F0 -F3  F2];...
           [ F2  F3  F0 -F1];...
           [ F3 -F2  F1  F0]];  % form approximation matrix explicitly once per patch

% evaluate basis at boundary nodes (to be multiplied by moments)
dx = rtau(1,:).*sbd.w; dy = rtau(2,:).*sbd.w; dz = rtau(3,:).*sbd.w;
[gradFbd,ijIdx] = evaltensorproductharmonicgrad_mex(r,order);   % harmonic gradient

% q^nm without kernel denominator
[qnm_0_1,qnm_0_2,qnm_0_3] = qnm_i(r,gradFbd,0);
[qnm_1_1,qnm_1_2,qnm_1_3] = qnm_i(r,gradFbd,1);
[qnm_2_1,qnm_2_2,qnm_2_3] = qnm_i(r,gradFbd,2);
[qnm_3_1,qnm_3_2,qnm_3_3] = qnm_i(r,gradFbd,3);

% omega^nm without kernel denominator
dr = [dx;dy;dz];
onm_0 = omeganm_i(r,dr,qnm_0_1,qnm_0_2,qnm_0_3);
onm_1 = omeganm_i(r,dr,qnm_1_1,qnm_1_2,qnm_1_3);
onm_2 = omeganm_i(r,dr,qnm_2_1,qnm_2_2,qnm_2_3);
onm_3 = omeganm_i(r,dr,qnm_3_1,qnm_3_2,qnm_3_3);

% moments (target dependent part)
if ~if_adapt
  M_all = momentsallplain_mex(r0,r,order,sbd.np,sbd.p);
  % M_all = momentsallplain(r0,r,order,sbd.np,sbd.p); % here to count flops
else
  M_all = momentsalladapt(r0,r,order,sbd.np,sbd.p);
end
% omega^nm with moments (target dependent part)
Omega_all = omegaall_mex(r0,M_all,order,reshape(onm_0,[],4),reshape(onm_1,[],4),reshape(onm_2,[],4),reshape(onm_3,[],4),h_dim,ijIdx);

% backwards-stable solve
warning('off','MATLAB:nearlySingularMatrix'); 
A = (Mmatrix'\Omega_all')'*Amu_c; % solve for special weights...

end

function M_all = momentsalladapt(r0,r,order,np,p)
% adaptive rule

rho = 4^(16/p);
[tj,wj,D0] = gauss(p);
w_bclag = bclag_interp_weights(tj);
Legmat = legendre.matrix(p); % for legendre coefficients
sqn_flag = false(np,numel(r0(1,:))); 
r_root = NaN(np,3,numel(r0(1,:))); t_root = NaN(np,numel(r0(1,:)));
rf_mex_flag = true; % rootfinder mex file flag
for ell = 1:np % loop over patch, find close targets, store
  idx_ell = (ell-1)*p + (1:p);
  r_ell = r(:,idx_ell);
  % find nearest pts on edge (designed to handle multiple targets)
  [sqn_idx_tmp,r_root_tmp,t_root_tmp] = line3_near_r0(tj, r_ell, r0, rho, Legmat,rf_mex_flag);
  % store roots info
  sqn_flag(ell,sqn_idx_tmp) = true;
  r_root(ell,:,sqn_idx_tmp) = r_root_tmp; t_root(ell,sqn_idx_tmp) = t_root_tmp;
end
[~,M_all] = momentsall_mex(r0,r,order,sqn_flag,r_root,t_root,tj,wj,D0,w_bclag,np,p);

end

function M_all = momentsallplain(r0,r,order,np,p)
% plain rule, matlab

M_all = zeros(numel(r(1,:)),2*order+2,numel(r0(1,:)));
for j=1:size(r0,2)
  % target
  r0_j = r0(:,j);
  % loop over edge patch
  for ell = 1:np
    idx_ell = (ell-1)*p + (1:p);
    r_ell = r(:,idx_ell);
    % plain smooth quadrature
    [~,mk_j] = moments(r0_j,r_ell,2*order+1);
    M_all(idx_ell,:,j) = mk_j;
  end
end

end

function s = get_high_order_quad(r,orderf)

% assume tensor product grid
order = sqrt(numel(r(1,:)));

% interpolation stuff
[x1 x2] = meshgrid(gauss(order)); 
uvs = [x1(:)';x2(:)'];
[xf,wf,D] = gauss(orderf);
[xf1 xf2] = meshgrid(xf);
uvsf = [xf1(:)';xf2(:)'];
[~,V,R] = interpmat_2d(uvsf,uvs);
wwf = wf(:)*wf(:)'; wwf = wwf(:)'; 

% input of grad SLP
sxf = ((V'\R')'*r')';
tp_ind = reshape(1:orderf^2,[orderf orderf]);
xts = nan(1,orderf^2); yts = xts; zts = xts; % partial t (small circle)
for j=1:orderf
  xk = sxf(1,tp_ind(:,j)); yk = sxf(2,tp_ind(:,j)); zk = sxf(3,tp_ind(:,j));
  xts(tp_ind(:,j)) = D*xk(:); yts(tp_ind(:,j)) = D*yk(:); zts(tp_ind(:,j)) = D*zk(:);
end
rts = [xts(:),yts(:),zts(:)]';
xps = nan(1,orderf^2); yps = xps; zps = xps; % partial p (big circle)
for j=1:orderf
  r = sxf(1,tp_ind(j,:)); y = sxf(2,tp_ind(j,:)); z = sxf(3,tp_ind(j,:));
  xps(tp_ind(j,:)) = D*r(:); yps(tp_ind(j,:)) = D*y(:); zps(tp_ind(j,:)) = D*z(:);
end
rps = [xps(:),yps(:),zps(:)]';
nx =  cross(rps, rts);       % outward normal
sp = sqrt(sum(nx.^2,1)); % speeds
nx = nx./sp;
w = sp.*wwf;

s.x = sxf;
s.nx = nx;
s.w = w;

end

function phi = outfun(a,b)
% define indicator function for torus
my_eps = 0.01; % if want to make this smaller, then need to refine reference soln
phi = @(x,y,z) ((z-a*cos(atan2(y,z))).^2 + (y-a*sin(atan2(y,z))).^2 + x.^2) > b^2+my_eps; 
end