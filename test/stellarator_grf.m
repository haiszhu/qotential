
if ~exist('rrqdir', 'var') || isempty(rrqdir)
  error('stellarator_grf:rrqdir', ...
        'Set rrqdir to the rrq checkout that holds test/stellarator (setup() and the .dat live there), then rerun.');
end
addpath(rrqdir);                                   % rrq-legacy root
addpath(fullfile(rrqdir, 'test', 'stellarator'));  % the .dat files
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'matlab'));  % qol_*
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
cd(fullfile(rrqdir, 'test', 'stellarator'));       % setup() and fopen are relative
setup()

isimd = 0;

ichart = 1;

MyOrder = [4 6 8 10 12 14 16];
MyMp = 12*[1 2 3 4 5 6 7 8];
MyNp = 36*[1 2 3 4 5 6 7 8];
MyTol = [1e-04 1e-06 1e-08 1e-09 1e-11 1e-14 1e-14];

JJ = [2 2 3 3 4 5 6];
MyErr = zeros(numel(JJ),1);
MyFMMtime = zeros(numel(JJ),1);
MyptswiseErr = [];

for idx = 1:numel(JJ)
  order = MyOrder(idx); 
  so.mp = MyMp(JJ(idx)); 
  so.np = MyNp(JJ(idx)); 

  [xg,wg,Dg] = gauss(order);
  [uvs,wts]  = get_vioreanu_nodes(order-1);
  [sx,snx,sw,rts,rps,npan] = stellarator_geo_mex(so.mp,so.np,order,xg,wg,Dg,uvs,wts);
  s.x = sx; s.nx = snx; s.w = sw; s.rts = rts; s.rps = rps;
  h_dim = order*(order+1)/2;

  pan = cell(1,npan);
  for k = 1:npan
    kk = (k-1)*h_dim + (1:h_dim);
    pan{k}.x   = sx(:,kk);
    pan{k}.nx  = snx(:,kk);
    pan{k}.w   = sw(kk);
    pan{k}.rts = rts(:,kk);
    pan{k}.rps = rps(:,kk);
  end

  [f, gradf]= torus_surf_den(); % for GRF
  ub = f(sx)'; ubn = sum(snx.*gradf(sx),1)';  % boundary u, and partial u partial normal
  sigma = sum(snx.*gradf(sx),1)'; tau = -f(sx)';     % densities for GRF (col vecs)
  
  t.x = s.x; % self eval
  t.nx = s.nx;
  fmm_eps = MyTol(idx);
  tic;
  sdlpval1 = Lap3dSLPfmm(t,s,ubn,fmm_eps) - Lap3dDLPfmm(t,s,ub,fmm_eps);
  fmmtime = toc;
  len = order*(order+1)/2;
  orderff = order+2;
  distff = 1.4;
  tmpidx = 1:numel(t.x(1,:));
  nquadadd = 2; 
  nquad = order + nquadadd;
  nterms = order;
  timeinfo = zeros(20,1);

  sbdnp = 3;  nbd = sbdnp*nquad;
  if ichart
    [xq,~,~]   = gauss(nquad);
    tpan = (0:sbdnp)'*2*pi/sbdnp;
    th1 = 2*pi/3;  th2 = 4*pi/3;  thi1 = 1/th1;
    uvbd = zeros(2,nbd+3);
    for e = 1:sbdnp
      for q = 1:nquad
        tt = tpan(e) + 0.5*(1+xq(q))*(tpan(e+1)-tpan(e));
        if tt < th1
          uvbd(:,(e-1)*nquad+q) = [0; 1 - thi1*tt];
        elseif tt < th2
          uvbd(:,(e-1)*nquad+q) = [thi1*(tt-th1); 0];
        else
          uvbd(:,(e-1)*nquad+q) = [1 - thi1*(tt-th2); thi1*(tt-th2)];
        end
      end
    end
    uvbd(:,nbd+1) = [0;0];  uvbd(:,nbd+2) = [1;0];  uvbd(:,nbd+3) = [0;1];
  end
  mp = so.mp; 
  np = so.np;
  hdim = order*(order+1)/2;
  npan = numel(pan);
  parameters = zeros(1,13);
  parameters(1:7) = [mp np nterms hdim nquad npan orderff];
  source = [sx;snx;sw;rts;rps]';
  data = [parameters;source];
  filename = sprintf('stellarator_source_mp%d_np%d_order%d.dat', mp, np, order);
  fid = fopen(filename, 'w');
  for i = 1:size(data, 1)
    fprintf(fid, '%.16e ', data(i, :));  
    fprintf(fid, '\n');  
  end
  fclose(fid);
  for k=1:numel(pan)
    
    qradii = 1.75*sqrt(sum(pan{k}.w)); % seems to be good
    qpoint = mean(pan{k}.x,2)';
    idxs = sum((t.x-qpoint').^2)<sum(qradii.^2);
    tc.x = t.x(:,idxs);
    idxs = tmpidx(idxs);
    ntc = length(tc.x(1,:));
  
    sk = pan{k};
    tx = tc.x;
    sxk = sk.x; swk = sk.w; snxk = sk.nx;
    rtsk = pan{k}.rts;
    rpsk = pan{k}.rps;
    m = size(tx,2);
    n = size(sxk,2);
    iside = 30; 
  
    S_ij = zeros(m,n);
    K_ij = zeros(m,n);
    IalphaAsvestas = zeros(m,1);
    
    order = 1/2*(-1 + sqrt(1+8*numel(sxk(1,:)))); % vioreanu-rokhlin nodes order
    nterms = order;
    hdim = order*(order+1)/2;
    len0 = 1;
    sbdnp = 3*len0;
    
    Omegas = zeros(4*hdim,m);
    
    distff = 1.4;
    IalphaAsvestas2 = zeros(size(IalphaAsvestas));
    if ichart
      xbuf = stellarator_geo_mex(so.mp,so.np,order,xg,wg,Dg,uvs,wts, k, uvbd);
    else
      xbuf = zeros(3,nbd+3);
    end
    [S_ij,K_ij,Omegas2,IalphaAsvestas2,timeinfo] = qol_rrq_mex(m,tx,iside,order,h_dim,nquad, ...
                             n,sxk,swk,snxk,rtsk,rpsk,orderff,distff,isimd, ...
                             ichart,xbuf(:,1:nbd),xbuf(:,nbd+1:nbd+3), ...
                             S_ij,K_ij,Omegas,IalphaAsvestas2,timeinfo);
  
    js = (k-1)*len+(1:len); 
    S_ij_naive = Lap3dSLPmat(tc,sk);
    K_ij_naive = Lap3dDLPmat(tc,sk);
    [~,selfFlag2] = ismember(idxs,js); % column: to mask source (was ismembc2, removed from MATLAB)
    j = nonzeros(selfFlag2)'; % column (idxs(something) shows up in js, its value is j)
    selfFlag = ismember(idxs,js); % row: to mask close targets (was ismembc, removed from MATLAB)
    tmp = 1:ntc; i = tmp(selfFlag); % row i(l), column j(l) is s.x(:,l)...
    selfIdx = i + (j-1)*ntc;
    S_ij_naive(selfIdx) = 0;
    K_ij_naive(selfIdx) = 0;
    sdlpval1(idxs) = sdlpval1(idxs) + ((S_ij-S_ij_naive)*ubn(js) - (K_ij-K_ij_naive)*ub(js));
  end
  gb1 = sdlpval1;
  err1 = abs(gb1)/max(abs(ubn));

  MyErr(idx) = max(err1(:));
  MyptswiseErr{idx} = err1;
  MyFMMtime(idx) = fmmtime;
  fprintf('case %d: order=%2d mp=%3d np=%3d  N=%d  GRF max rel err = %.3e  (fmm %.2f s)\n', ...
          idx, order, so.mp, so.np, numel(sx(1,:)), MyErr(idx), fmmtime);
  save('stellarator.mat',"MyErr","MyptswiseErr","MyFMMtime");
end

function [sl,su,s_sub] = get_high_order_quad(pan,order,qntype,side)
[x,w,D]=gauss(order); [x1 x2] = meshgrid(x); ww = w(:)*w(:)'; 
uvs = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
wts = ww(:)';
sx = pan.x;
if side == 'i'
  sx1 = reshape(sx(1,:),order,order);
  sx2 = reshape(sx(2,:),order,order);
  sx3 = reshape(sx(3,:),order,order);
  sx1 = sx1(:,end:-1:1);
  sx2 = sx2(:,end:-1:1);
  sx3 = sx3(:,end:-1:1);
  sx = [sx1(:) sx2(:) sx3(:)]';
end
tp_ind = reshape(1:order^2,[order order]);
xts = nan(1,order^2); yts = xts; zts = xts; % partial t (small circle)
for j=1:order
  xk = sx(1,tp_ind(:,j)); yk = sx(2,tp_ind(:,j)); zk = sx(3,tp_ind(:,j));
  xts(tp_ind(:,j)) = D*xk(:); yts(tp_ind(:,j)) = D*yk(:); zts(tp_ind(:,j)) = D*zk(:);
end
rts = [xts(:),yts(:),zts(:)]';
xps = nan(1,order^2); yps = xps; zps = xps; % partial p (big circle)
for j=1:order
  x = sx(1,tp_ind(j,:)); y = sx(2,tp_ind(j,:)); z = sx(3,tp_ind(j,:));
  xps(tp_ind(j,:)) = D*x(:); yps(tp_ind(j,:)) = D*y(:); zps(tp_ind(j,:)) = D*z(:);
end
rps = [xps(:),yps(:),zps(:)]';
snx =  cross(rps, rts); % outward normal
ssp = sqrt(sum(snx.^2,1)); % speeds
snx = snx./ssp;  % surface normal
sw = ssp.*wts;  % quadrature weights

sl.x = sx; sl.nx = snx; sl.w = sw;

su = [];
s_sub = [];
[xc yc] = meshgrid([-1/2 1/2]); factor = 1/2;
for k=1:numel(xc)
  xck = xc(k); yck = yc(k);
  uvsk = factor*[x1(:)';x2(:)']+[xck;yck];
  [Lk,Vk,Rk] = interpmat_2d(uvsk,uvs);
  sxk = (Rk*(Vk\sx'))';
  xtsk = nan(1,order^2); ytsk = xtsk; ztsk = xtsk; % partial t (small circle)
  for j=1:order
    xk = sxk(1,tp_ind(:,j)); yk = sxk(2,tp_ind(:,j)); zk = sxk(3,tp_ind(:,j));
    xts(tp_ind(:,j)) = D*xk(:); yts(tp_ind(:,j)) = D*yk(:); zts(tp_ind(:,j)) = D*zk(:);
  end
  rtsk = [xts(:),yts(:),zts(:)]';
  xpsk = nan(1,order^2); ypsk = xpsk; zpsk = xpsk; % partial p (big circle)
  for j=1:order
    x = sxk(1,tp_ind(j,:)); y = sxk(2,tp_ind(j,:)); z = sxk(3,tp_ind(j,:));
    xps(tp_ind(j,:)) = D*x(:); yps(tp_ind(j,:)) = D*y(:); zps(tp_ind(j,:)) = D*z(:);
  end
  rpsk = [xps(:),yps(:),zps(:)]';
  snxk =  cross(rpsk, rtsk); % outward normal
  sspk = sqrt(sum(snxk.^2,1)); % speeds
  snxk = snxk./sspk;  % surface normal
  swk = sspk.*wts;  % quadrature weights
  s_sub{k}.x = rpsk; s_sub{k}.nx = snxk; s_sub{k}.w = swk; 
end
end

function [f,gradf] = torus_surf_den()
lam = 0.5; f = @(x) exp(lam*x(1,:)).*cos(lam*x(2,:));
gradf = @(x) [lam*f(x); -lam*exp(lam*x(1,:)).*sin(lam*x(2,:)); 0*f(x)];
end

function pan3 = stellarator_pan_init(pan0,order,curv_tol)

pan = pan0; 

pan2 = []; addind = 0;
[x0,~,~] = gauss(order); [x1 x2] = meshgrid(x0); xx = [x1(:)';x2(:)']; 
for k=1:numel(pan)
  dist_t = norm(pan{k}.x(:,1) - pan{k}.x(:,end-order+1)); % along t-direction, if I recall correctly...
  dist_p = norm(pan{k}.x(:,1) - pan{k}.x(:,order));       % along p-direction
  if dist_t/dist_p > 1.65
    orig = pan0{k}.para(1);
    j = pan0{k}.para(2);
    i = pan0{k}.para(3);
    tpansiz1 = pan0{k}.para(4);
    tpansiz2 = pan0{k}.para(5);
    pank1.chart = @(t) stellaratorparam(((t(1,:)/2-1/2)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2); 
    pank2.chart = @(t) stellaratorparam(((t(1,:)/2+1/2)+orig-1+2*j)*tpansiz1,(t(2,:)+orig-1+2*i)*tpansiz2); 
    pan2{k+addind} = pank1;
    pan2{k+addind+1} = pank2;
    addind = addind+1;
  elseif dist_p/dist_t > 1.65
    orig = pan0{k}.para(1);
    j = pan0{k}.para(2);
    i = pan0{k}.para(3);
    tpansiz1 = pan0{k}.para(4);
    tpansiz2 = pan0{k}.para(5);
    pank1.chart = @(t) stellaratorparam((t(1,:)+orig-1+2*j)*tpansiz1,((t(2,:)/2-1/2)+orig-1+2*i)*tpansiz2); 
    pank2.chart = @(t) stellaratorparam((t(1,:)+orig-1+2*j)*tpansiz1,((t(2,:)/2+1/2)+orig-1+2*i)*tpansiz2);
    pan2{k+addind} = pank1;
    pan2{k+addind+1} = pank2;
    addind = addind+1;
  else
    pan2{k+addind} = pan{k};
  end
end
for k = 1:numel(pan2)
  sxk = pan2{k}.chart(xx);
  pan2{k} = get_gl_quadr(sxk,order);
end
pan3 = pan2;

end
