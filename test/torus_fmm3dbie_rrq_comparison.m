% Compare RRQ and FMM3DBIE near quadrature on the same torus problem.

clear

rrq_total_tic = tic;
fprintf('\n============================================================\n');
fprintf('Torus near-quadrature comparison\n');
fprintf('============================================================\n');

profile clear
profile on

qroot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(qroot,'matlab'))
addpath(fullfile(qroot,'utils'))
addpath(fullfile(qroot,'external','LineQuaaadrature','utils'))
addpath(fullfile(qroot,'external','LineQuaaadrature','matlab'))
addpath(fullfile(qroot,'external','QuatApproximation','matlab'))
addpath(fullfile(qroot,'external','kdtree','toolbox'))

order = 8;
mp = 2*6;
np = 2*10;
orig = 0;
radii_rrq = [1/2; 1; 0];
scales_rrq = [1; 1; 1];
nosc_rrq = 0;
hdim = order*(order+1)/2;
ntri = 2*mp*np;
nsrc = ntri*hdim;
geometry_tic = tic;
[sx,snx,sw,rts,rps,ier] = lqtm_create_torus_tri_mesh_mex( mp,np,order,radii_rrq,scales_rrq,nosc_rrq,orig,ntri, ...
                                                          zeros(3,nsrc),zeros(3,nsrc),zeros(nsrc,1), ...
                                                          zeros(3,nsrc),zeros(3,nsrc),0);
assert(ier == 0, 'lqtm_create_torus_tri_mesh_mex failed with ier=%d',ier);
sw = sw(:).';
pan = cell(ntri,1);
for k = 1:ntri
    idx = (k-1)*hdim+(1:hdim);
    pan{k}.x = sx(:,idx);
    pan{k}.nx = snx(:,idx);
    pan{k}.w = sw(idx);
    pan{k}.rts = rts(:,idx);
    pan{k}.rps = rps(:,idx);
end
s.x = sx; s.nx = snx; s.w = sw; s.rts = rts; s.rps = rps;
geometry_time = toc(geometry_tic);

fprintf('\n[RRQ]\n');
fprintf('  geometry   : order=%d, panels=%d, nodes=%d, area=%.15e, time=%.3f s\n',order,numel(pan),size(s.x,2),sum(s.w),geometry_time);

% targets
t.x = s.x; % self eval
t.nx = s.nx;

naive_matrix_tic = tic;
AmatS = zeros(numel(s.x(1,:)),numel(s.x(1,:)));
AmatS(:,1:end/8) = Lap3dSLPmat(t,struct('x',s.x(:,1:end/8),'nx',s.nx(:,1:end/8),'w',s.w(1:end/8)));
AmatS(:,end/8+1:end/4) = Lap3dSLPmat(t,struct('x',s.x(:,end/8+1:end/4),'nx',s.nx(:,end/8+1:end/4),'w',s.w(end/8+1:end/4)));
AmatS(:,end/4+1:3*end/8) = Lap3dSLPmat(t,struct('x',s.x(:,end/4+1:3*end/8),'nx',s.nx(:,end/4+1:3*end/8),'w',s.w(end/4+1:3*end/8)));
AmatS(:,3*end/8+1:end/2) = Lap3dSLPmat(t,struct('x',s.x(:,3*end/8+1:end/2),'nx',s.nx(:,3*end/8+1:end/2),'w',s.w(3*end/8+1:end/2)));
AmatS(:,end/2+1:5*end/8) = Lap3dSLPmat(t,struct('x',s.x(:,end/2+1:5*end/8),'nx',s.nx(:,end/2+1:5*end/8),'w',s.w(end/2+1:5*end/8)));
AmatS(:,5*end/8+1:3*end/4) = Lap3dSLPmat(t,struct('x',s.x(:,5*end/8+1:3*end/4),'nx',s.nx(:,5*end/8+1:3*end/4),'w',s.w(5*end/8+1:3*end/4)));
AmatS(:,3*end/4+1:7*end/8) = Lap3dSLPmat(t,struct('x',s.x(:,3*end/4+1:7*end/8),'nx',s.nx(:,3*end/4+1:7*end/8),'w',s.w(3*end/4+1:7*end/8)));
AmatS(:,7*end/8+1:end) = Lap3dSLPmat(t,struct('x',s.x(:,7*end/8+1:end),'nx',s.nx(:,7*end/8+1:end),'w',s.w(7*end/8+1:end)));
AmatD = zeros(numel(s.x(1,:)),numel(s.x(1,:)));
AmatD(:,1:end/8) = Lap3dDLPmat(t,struct('x',s.x(:,1:end/8),'nx',s.nx(:,1:end/8),'w',s.w(1:end/8)));
AmatD(:,end/8+1:end/4) = Lap3dDLPmat(t,struct('x',s.x(:,end/8+1:end/4),'nx',s.nx(:,end/8+1:end/4),'w',s.w(end/8+1:end/4)));
AmatD(:,end/4+1:3*end/8) = Lap3dDLPmat(t,struct('x',s.x(:,end/4+1:3*end/8),'nx',s.nx(:,end/4+1:3*end/8),'w',s.w(end/4+1:3*end/8)));
AmatD(:,3*end/8+1:end/2) = Lap3dDLPmat(t,struct('x',s.x(:,3*end/8+1:end/2),'nx',s.nx(:,3*end/8+1:end/2),'w',s.w(3*end/8+1:end/2)));
AmatD(:,end/2+1:5*end/8) = Lap3dDLPmat(t,struct('x',s.x(:,end/2+1:5*end/8),'nx',s.nx(:,end/2+1:5*end/8),'w',s.w(end/2+1:5*end/8)));
AmatD(:,5*end/8+1:3*end/4) = Lap3dDLPmat(t,struct('x',s.x(:,5*end/8+1:3*end/4),'nx',s.nx(:,5*end/8+1:3*end/4),'w',s.w(5*end/8+1:3*end/4)));
AmatD(:,3*end/4+1:7*end/8) = Lap3dDLPmat(t,struct('x',s.x(:,3*end/4+1:7*end/8),'nx',s.nx(:,3*end/4+1:7*end/8),'w',s.w(3*end/4+1:7*end/8)));
AmatD(:,7*end/8+1:end) = Lap3dDLPmat(t,struct('x',s.x(:,7*end/8+1:end),'nx',s.nx(:,7*end/8+1:end),'w',s.w(7*end/8+1:end)));
naive_matrix_time = toc(naive_matrix_tic);
naive_matrix_gib = 2*8*size(s.x,2)^2/1024^3;

fprintf('  dense      : storage=%.3f GiB, time=%.3f s\n',naive_matrix_gib,naive_matrix_time);

% compute close targs
kd = KDTree(t.x'); 
% loop over patch and compute close interaction nodes
tcx = zeros(size(t.x));
tcxi = ones(numel(pan)+1,1);
tcxrow = zeros(size(t.x,2),1);
close_search_tic = tic;
for k=1:numel(pan)
  % targets close to this patch
  qradii = 1.75*sqrt(sum(pan{k}.w)); % seems to be good
  qpoint = mean(pan{k}.x,2)';
  idxcin = kd.ball( qpoint, qradii); 
  tcjin.x = t.x(:,idxcin);
  %
  tcxi(k+1) = tcxi(k) + size(tcjin.x,2);
  tcx(:,tcxi(k):(tcxi(k+1)-1)) = tcjin.x;
  tcxrow(tcxi(k):(tcxi(k+1)-1)) = idxcin;
end
tcx = tcx(:,1:tcxi(end)-1);
tcxrow = tcxrow(1:tcxi(end)-1);
close_search_time = toc(close_search_tic);
close_counts = diff(tcxi);

fprintf('  near search: pairs=%d, targets/panel=%d/%.2f/%d, time=%.3f s\n',size(tcx,2),min(close_counts),mean(close_counts),max(close_counts),close_search_time);

timeinfo = zeros(20,1);
nterms = order;
hdim = order*(order+1)/2;
npan = numel(pan);
nquad = nterms + 2;
orderff = nterms + 2;

tmpidxin = 1:numel(t.x(1,:));
korder = order-1;
kpols = order*(order+1)/2;
vmatr = zeros(kpols,kpols);
umatr = zeros(kpols,kpols);
[umatr,vmatr] = qakg_koorn_vals2coefs_coefs2vals_mex(korder,kpols,umatr,vmatr);

%
npan = numel(pan);
nsrc = size(sx,2);
ntcx = size(tcx,2);

% we know sparse matrix nnz info
nnz = (tcxi(end)-1)*hdim; 
sparsei = zeros(nnz,1);
sparsej = zeros(nnz,1);
sparsevs = zeros(nnz,1);
sparsevd = zeros(nnz,1);

%
distff = 1.4;
iside = 30; 
isimd = 0;

%% rrq version
rrq_tic = tic;
[sparsei, sparsej, sparsevs, sparsevd] = qol_getnearquad_lap_rrq_mex(npan, nterms, hdim, nquad, ...
                            nsrc, sx, snx, sw, rts, rps, ... % source info
                            ntcx, tcx, tcxrow, tcxi, ... % target info
                            orderff, distff, iside, isimd, ... % parameters
                            nnz, sparsei, sparsej, sparsevs, sparsevd, timeinfo);
rrq_time = toc(rrq_tic);
rrq_pps = nsrc/rrq_time;

fprintf('  correction : nquad=%d, orderff=%d, distff=%.2f, entries=%d, time=%.3f s, PPS=%.3e\n',nquad,orderff,distff,nnz,rrq_time,rrq_pps);

sparse_assembly_tic = tic;
AtcxmatS = sparse(sparsei,sparsej,sparsevs,numel(s.x(1,:)),numel(s.x(1,:)));
AtcxmatD = sparse(sparsei,sparsej,sparsevd,numel(s.x(1,:)),numel(s.x(1,:)));
sparse_assembly_time = toc(sparse_assembly_tic);

if 0

  %
  formulation = 'SLP + DLP';
  AmatS(diagind(AmatS)) = 0;
  AmatD(diagind(AmatD)) = 0;
  Amat = AmatS + AmatD + AtcxmatS + AtcxmatD;
  
  y_source.x = [1.0;0;-0.02]; y_source.w = 1; % source inside torus
  rng(4)
  pt_source = rand(1,1); pt_source = 10*pt_source/norm(pt_source);
  A_rhs = Lap3dSLPmat(s,y_source);
  rhs = A_rhs*pt_source;
  
  % solve
  linear_solve_tic = tic;
  tau = Amat\rhs;
  linear_solve_time = toc(linear_solve_tic);
  
  % target
  tout.x = [1;1;0.6];
  target_eval_tic = tic;
  Atmat = Lap3dSLPmat(tout,s) + Lap3dDLPmat(tout,s);
  val = Atmat*tau;
  target_eval_time = toc(target_eval_tic);

else
  
  %
  formulation = 'SLP';
  AmatS(diagind(AmatS)) = 0;
  Amat = AmatS + AtcxmatS;
  
  y_source.x = [1.0;0;-0.02]; y_source.w = 1; % source inside torus
  rng(4)
  pt_source = rand(1,1); pt_source = 10*pt_source/norm(pt_source);
  A_rhs = Lap3dSLPmat(s,y_source);
  rhs = A_rhs*pt_source;
  
  % solve
  linear_solve_tic = tic;
  tau = Amat\rhs;
  linear_solve_time = toc(linear_solve_tic);
  
  % target
  tout.x = [1;1;0.6];
  target_eval_tic = tic;
  Atmat = Lap3dSLPmat(tout,s);
  val = Atmat*tau;
  target_eval_time = toc(target_eval_tic);

end

%
reference_eval_tic = tic;
Ahom = Lap3dSLPmat(tout,y_source);
refval = Ahom*pt_source;
reference_eval_time = toc(reference_eval_tic);

abs_error = abs(val-refval);
rel_error = abs_error/max(abs(refval),realmin);
rrq_total_time = toc(rrq_total_tic);

fprintf('  solve      : formulation=%s, time=%.3f s, relative error=%.3e, total=%.3f s\n',formulation,linear_solve_time,rel_error,rrq_total_time);

% keyboard

%% fmm3dbie version (modern vs legacy)

% if_modern = true; % something is still not quite right the way I used modern inferface, more than just change of src ordering
if_modern = false;

% Current official fmm3dbie MATLAB interface.
run(fullfile(getenv('HOME'),'git','fmm3dbie','matlab','startup.m'));
% Match the RRQ order-6 discretization:
% fmm3dbie's triangular patch order is one less than RRQ's node count order.
norder = order - 1;
% Historical get_wtorus_geom used:
%   radii = [minor radius; major radius; wave amplitude]
radii = [1/2; 1; 0];
scales = [1; 1; 1];
nosc = 0; % w
nu = 2*10;
nv = 2*6;

fmm3dbie_total_tic = tic;
fmm3dbie_geometry_tic = tic;
if if_modern
  % geometries.startorus instead expects:
  %   [major radius; minor radius; wave amplitude]
  S_fmm = geometries.startorus( radii([2,1,3]), nosc, scales, [nu,nv], norder, 1);
  [srcvals, srccoefs, norders, ixyzs, iptype, wts] = extract_arrays(S_fmm);
  npatches = S_fmm.npatches;
  npts = S_fmm.npts;
else
  npatches = 2*nu*nv;
  % call
  [norders,ixyzs,iptype,srcvals,srccoefs,wts] = get_wtorus_geom_mex(radii,scales,nosc,nu,nv,norder);
  npatches = numel(norders);
  npts = size(srcvals,2);
end
fmm3dbie_geometry_time = toc(fmm3dbie_geometry_tic);
assert(npatches == 2*nu*nv);
assert(npts == npatches*order*(order+1)/2);

% Plot geometry.
figure(1), clf
plot3(srcvals(1,:), srcvals(2,:), srcvals(3,:), '.');
axis equal
hold on
area_exact = 4*pi^2*radii(1)*radii(2);
fprintf('\n[FMM3DBIE]\n');
fprintf('  geometry   : order=%d, panels=%d, nodes=%d, area error=%.3e, time=%.3f s\n',order,npatches,npts,abs(sum(wts)-area_exact),fmm3dbie_geometry_time);


% Laplace combined-field parameters:
% alpha = SLP strength; beta = DLP strength.
alpha = 1.0;
if 0
  %
  formulation = 'SLP + DLP';
  beta = 1.0;
else
  %
  formulation = 'SLP';
  beta = 0.0;
end
dpars = [alpha; beta];

eps_fmm3dbie = 0.51e-8;

% Build the complete fmm3dbie near-quadrature correction.
fmm3dbie_near_tic = tic;
if if_modern
  opts_quad = struct('format', 'rsc');
  Q_fmm = lap3d.dirichlet.get_quadrature_correction( ...
      S_fmm, eps_fmm3dbie, dpars, S_fmm, opts_quad);
  % Retain the old variable names where useful for inspection.
  row_ptr = Q_fmm.row_ptr;
  col_ind = Q_fmm.col_ind;
  iquad = Q_fmm.iquad;
  wnear = Q_fmm.wnear;
  nquad_fmm3dbie = Q_fmm.nquad;
else
  npts = numel(srcvals(1,:)); % srcvals are targets (when set up BIE boundary condition)
  % breakdown of lap_comb_dir_solver
  pi = atan(1)*4;
  % c
  % c
  % c        setup targets as on surface discretization points
  % c 
  ndtarg = 3;
  ntarg = npts;
  targs=zeros(ndtarg,npts); uvs_targ=zeros(2,ntarg); ipatch_id=zeros(ntarg,1);
  for i=1:ntarg
    targs(1,i) = srcvals(1,i);
    targs(2,i) = srcvals(2,i);
    targs(3,i) = srcvals(3,i);
    ipatch_id(i) = -1;
    uvs_targ(1,i) = 0;
    uvs_targ(2,i) = 0;
  end
  npatchesp1 = npatches + 1;
  mex_id_ = ['get_patch_id_uvs(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io double[xx])'];
  % [ipatch_id,uvs_targ] = get_patch_id_uvs_mex(npatches,norders,ixyzs,iptype,npts,ipatch_id,uvs_targ);
  [ipatch_id,uvs_targ] = fmm3dbie_routs(mex_id_,npatches,norders,ixyzs,iptype,npts,ipatch_id,uvs_targ,1,npatches,npatchesp1,npatches,1,npts,2,npts);
  iptype_avg = floor(sum(iptype)/(npatches+0.0d0));
  norder_avg = floor(sum(norders)/(npatches+0.0d0));
  rfac = 0.0; rfac0 = 0.0;
  mex_id_ = 'get_rfacs(c i int64_t[x], c i int64_t[x], c io double[x], c io double[x])';
  % [rfac,rfac0] = get_rfacs_mex(norder_avg,iptype_avg); % estimates near field parameters
  [rfac,rfac0] = fmm3dbie_routs(mex_id_,norder_avg,iptype_avg,rfac,rfac0,1,1,1,1);
  cms=zeros(3,npatches); rads=zeros(npatches,1); rad_near=zeros(npatches,1);
  
  [cms,rads] = get_centroid_rads_mex(npatches,norders,ixyzs,iptype,npts,srccoefs); % compute the centroid and bounding sphere radii for a collection of patches (calls get_centroid_rads_tri)
  
  for i=1:npatches
    rad_near(i) = rads(i)*rfac;
  end
  
  nnz = 0;
  mex_id_ = 'findnearmem(c i double[xx], c i int64_t[x], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c io int64_t[x])';
  % nnz = findnearmem_mex(cms,npatches,rad_near,targs,npts); % find near quadrature correction interactions
  [nnz] = fmm3dbie_routs(mex_id_,cms,npatches,rad_near,ndtarg,targs,ntarg,nnz,3,npatches,1,npatches,1,ndtarg,ntarg,1,1);
  
  row_ptr=zeros(npts+1,1); col_ind=zeros(nnz,1);
  mex_id_ = 'findnear(c i double[xx], c i int64_t[x], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c io int64_t[x], c io int64_t[x])';
  % [row_ptr,col_ind] = findnear_mex(cms,npatches,rad_near,targs,npts,nnz); %
  [row_ptr,col_ind] = fmm3dbie_routs(mex_id_,cms,npatches,rad_near,ndtarg,targs,ntarg,row_ptr,col_ind,3,npatches,1,npatches,1,ndtarg,ntarg,1,ntarg+1,nnz);
  
  %
  iquad = zeros(nnz+1,1);
  mex_id_ = 'get_iquad_rsc(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x])';
  % iquad = get_iquad_rsc_mex(npatches,ixyzs,npts,nnz,row_ptr,col_ind);
  [iquad] = fmm3dbie_routs(mex_id_,npatches,ixyzs,npts,nnz,row_ptr,col_ind,iquad,1,npatchesp1,1,1,npts+1,nnz,nnz+1);
  
  % %
  % ikerorder = -1;
  % if(abs(dpars(2))>1.0d-16), ikerorder = 0; end
  % novers=zeros(npatches,1);
  % ixyzso=zeros(npatches+1,1);
  % %
  % ztmp = 0;
  % mex_id_ = 'get_far_order(c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i dcomplex[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c io int64_t[x], c io int64_t[x])';
  % % [novers, ixyzso] = get_far_order_mex(eps,npatches,norders,ixyzs,iptype,cms,rads,npts,srccoefs,npts,targs,ikerorder,nnz,row_ptr,col_ind,rfac); % estimate oversampling for far-field, and oversample geometry
  % [novers,ixyzso] = fmm3dbie_routs(mex_id_,eps_fmm3dbie,npatches,norders,ixyzs,iptype,cms,rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,ztmp,nnz,row_ptr,col_ind,rfac,novers,ixyzso,1,1,npatches,npatchesp1,npatches,3,npatches,npatches,1,9,npts,1,1,ndtarg,ntarg,1,1,1,ntarg+1,nnz,1,npatches,npatchesp1);
  % npts_over = ixyzso(npatches+1)-1;
  % srcover=zeros(12,npts_over); wover=zeros(npts_over,1);
  % srcover = oversample_geom_mex(npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,novers,ixyzso,npts_over); % This subroutine oversamples geometry information given expansion coeffs of geometry info on patches
  % wover = get_qwts_mex(npatches,novers,ixyzso,iptype,npts_over,srcover); % smooth quadrature weights at the discretization nodes
  
  % c
  % c   compute near quadrature correction
  % c
  nquad = iquad(nnz+1)-1;
  % print *, "nquad=",nquad
  wnear = zeros(nquad,1);
        
  iquadtype = 1;
  % print *, "starting to generate near quadrature"
  mex_id_ = 'getnearquad_lap_comb_dir_eval(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c io double[x])';
  % wnear = getnearquad_lap_comb_dir_mex(npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,eps,dpars,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear);
  [wnear] = fmm3dbie_routs(mex_id_,npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,ipatch_id,uvs_targ,eps_fmm3dbie,dpars,iquadtype,nnz,row_ptr,col_ind,iquad,rfac0,nquad,wnear,1,npatches,npatchesp1,npatches,1,9,npts,12,npts,1,1,ndtarg,ntarg,ntarg,2,ntarg,1,2,1,1,ntarg+1,nnz,nnz+1,1,1,nquad);
  nquad_fmm3dbie = numel(wnear);
end
fmm3dbie_near_time = toc(fmm3dbie_near_tic);
nnz_fmm3dbie = numel(col_ind);
fmm3dbie_pps = npts/fmm3dbie_near_time;

fprintf('  correction : eps=%.2e, nnz=%d, nquad=%d, time=%.3f s, PPS=%.3e\n',eps_fmm3dbie,nnz_fmm3dbie,nquad_fmm3dbie,fmm3dbie_near_time,fmm3dbie_pps);

% max number of iterations
numit = 200;

if if_modern
  s2.x = srcvals(1:3,:);
  s2.nx = srcvals(10:12,:);
  s2.w = wts(:).';
  t2.x = s2.x;

  Amat2 = Lap3dSLPmat(t2,s2);
  Amat2(diagind(Amat2)) = 0;

  rhs2 = Lap3dSLPmat(s2,y_source)*pt_source;
  Atmat2 = Lap3dSLPmat(tout,s2);
else
  Amat2 = AmatS;
  rhs2 = rhs;
  Atmat2 = Atmat;
end

for itarg = 1:npts
  for j = row_ptr(itarg):row_ptr(itarg+1)-1
    ipatch = col_ind(j);
    source_ids = ixyzs(ipatch):ixyzs(ipatch+1)-1;
    weight_ids = iquad(j):iquad(j+1)-1;
    Amat2(itarg,source_ids) = wnear(weight_ids).';
  end
end
fmm3dbie_solve_tic = tic;
tau2 = Amat2\rhs2;
fmm3dbie_solve_time = toc(fmm3dbie_solve_tic);

val2 = Atmat2*tau2;
abs_error2 = abs(val2-refval);
rel_error2 = abs_error2/max(abs(refval),realmin);
fmm3dbie_total_time = toc(fmm3dbie_total_tic);

fprintf('  solve      : formulation=%s, time=%.3f s, relative error=%.3e, total=%.3f s\n',formulation,fmm3dbie_solve_time,rel_error2,fmm3dbie_total_time);

profile viewer


% keyboard


function [sl,su,s_sub] = get_high_order_quad(pan,order,qntype)
 [x,w]=gauss(order); [x1 x2] = meshgrid(x); ww = w(:)*w(:)'; 
 uvs = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
 wts = ww(:)';
 % 
 [rs, nrs, sps] = pan.chart(uvs); % these should be interpolated, update 'panel_smooth_quadk.m'
 ws = sps.*wts*pan.spfac;
 sl.x = rs; sl.nx = nrs; sl.w = ws; 
 sl.chart = pan.chart;
%  sl.bdchart = pan.bdchart; 
 su = [];
 %
 s_sub = [];
 [xc yc] = meshgrid([-1/2 1/2]); factor = 1/2;
 for k=1:numel(xc)
   xck = xc(k); yck = yc(k);
   uvsk = factor*[x1(:)';x2(:)']+[xck;yck];
   % 
   [rs, nrs, sps] = pan.chart(uvsk); % these should be interpolated, update 'panel_smooth_quadk.m'
   ws = sps.*wts*pan.spfac/4;
   s_sub{k}.x = rs; s_sub{k}.nx = nrs; s_sub{k}.w = ws; 
   s_sub{k}.spfac = pan.spfac/4;
   s_sub{k}.chart =@(t) pan.chart(factor*t+[xck;yck]);
 end
end

function phi = inoutfun(so,side)
% define indicator function for torus
my_eps = 0.125;
my_eps = 0.05;
if side == 'e'
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) > so.b^2+my_eps; 
else
  phi = @(x,y,z) ((x-so.a*cos(atan2(y,x))).^2 + (y-so.a*sin(atan2(y,x))).^2 + z.^2) < so.b^2-my_eps;
end
end
