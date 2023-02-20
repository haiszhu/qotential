function [s, N, np] = quadr_3dline(s, N, qntype, xwD)  % set up quadrature on a closed segment
% similar to QUADR, more later
% Hai 08/31/22, add support to derive panel boundary map from provided
% vioreanu nodes, in case s.Z is not available. 09/14/22


if nargin==0, test_quadr_3dline; return; end
if ~isfield(s,'p'), s.p=16; end, p = s.p;
if nargin < 4, 
    if qntype=='G', [x, w, D] = gauss(p); else, [x, w, D] = cheby(p); end, 
else
    x = xwD.x; w = xwD.w; D = xwD.D; 
end

if ~isfield(s,'tpan')
    np = ceil(N/p); N = p*np;      % np = # panels
    s.tlo = (0:np-1)/np*2*pi; % panel start params
    s.thi = (1:np)/np*2*pi; % panel end params
    s.np = np;
else
    np = numel(s.tpan)-1; N = p*np;
    s.tlo = s.tpan(1:end-1); 
    s.thi = s.tpan(2:end);
    s.np = np;
end

if isfield(s,'Z')
  s.xlo = s.Z(s.tlo);
  s.xhi = s.Z(s.thi);
else % only support on-surface vioreanu nodes & tensor product G-L nodes
  x_uvs = s.xvr; % will change field name: vr was short for vioreanu-rokhlin
  order = 1/2*(-1 + sqrt(1+8*numel(x_uvs(1,:)))); % vioreanu-rokhlin nodes
  ordert = sqrt(numel(x_uvs(1,:)));
  if isfield(s,'x_uvs_qntype')
    x_uvs_qntype = s.x_uvs_qntype; 
  else
    x_uvs_qntype = 'V'; % to be compatible with old code, Stoks etc. which only support vioreanu-rokhlin
  end
  if x_uvs_qntype == 'V'
    n_idx = 0:order-1; k_idx = 0:order-1;
    [ntmp,ktmp] = meshgrid(n_idx,k_idx); uptri_idx = logical(triu(ones(order)));
    Kvr = ktmp(uptri_idx); Nvr = ntmp(uptri_idx);
    polsKN = [Kvr,Nvr];
    if isfield(s,'umatr')
      umatr = s.umatr;
    else
      [uvs,wts]=get_vioreanu_nodes(order-1); 
      [umatr,~]=koorn_vals2coefs_coefs2vals(order-1,order*(order+1)/2,uvs,polsKN);
    end
    coefs_x = umatr*x_uvs(1,:)'; coefs_y = umatr*x_uvs(2,:)'; coefs_z = umatr*x_uvs(3,:)';
    striZ = @(t) (mod(t,2*pi)<pi/2).*(0+1i*(1-2/pi*mod(t,2*pi)))...
              +((mod(t,2*pi)>=pi/2)&(mod(t,2*pi)<pi)).*(2/pi*(mod(t,2*pi)-pi/2)+0*1i)...
              +(mod(t,2*pi)>=pi).*(1-1/pi*(mod(t,2*pi)-pi)+1i*(1/pi*(mod(t,2*pi)-pi)));
    coefs_xyz = [coefs_x,coefs_y,coefs_z];
    if isfield(s,'xlo_mat')
      s.xlo = (s.xlo_mat*coefs_xyz)';
    else
      s.xlo = (koorn_coefs2vals_vec(order-1,order*(order+1)/2,[real(striZ(s.tlo));imag(striZ(s.tlo))],polsKN)*coefs_xyz)';
    end
    if isfield(s,'xhi_mat')
      s.xhi = (s.xhi_mat*coefs_xyz)';
    else
      s.xhi = (koorn_coefs2vals_vec(order-1,order*(order+1)/2,[real(striZ(s.thi));imag(striZ(s.thi))],polsKN)*coefs_xyz)';
    end
    s.trichart = @(t) (koorn_coefs2vals_vec(order-1,order*(order+1)/2,t,polsKN)*coefs_xyz)';
  elseif x_uvs_qntype == 'T' % tensor product grid
    xtmp = gauss(ordert);   % Gauss-Legendre nodes & weights on [-1,1]
    [x1 x2] = meshgrid(xtmp); xx = [x1(:)';x2(:)'];  % 2*p^2 parameter col vecs in R2
    rboxZ = @(t) (mod(t,2*pi)<pi/2).*(-1+4/pi*mod(t,2*pi)-1i)...
            +((mod(t,2*pi)>=pi/2)&(mod(t,2*pi)<pi)).*(1+1i*(-1+4/pi*(mod(t,2*pi)-pi/2)))...
            +((mod(t,2*pi)>=pi)&(mod(t,2*pi)<3*pi/2)).*(1-4/pi*(mod(t,2*pi)-pi)+1i)...
            +(mod(t,2*pi)>=3*pi/2).*(-1+1i*(1-4/pi*(mod(t,2*pi)-3*pi/2)));
    if isfield(s,'xlo_mat')
      s.xlo = (s.xlo_mat*x_uvs')';
    else
      s.xlo = interpval([real(rboxZ(s.tlo));imag(rboxZ(s.tlo))],x_uvs,xtmp);
    end
    if isfield(s,'xhi_mat')
      s.xhi = (s.xhi_mat*x_uvs')';
    else
      s.xlo = interpval([real(rboxZ(s.thi));imag(rboxZ(s.thi))],x_uvs,xtmp);
    end
    s.quadchart = @(t) (interpmat_2d(t,xx)*x_uvs')';
    s.quadchart = @(t) interpval(t,x_uvs,xtmp);
  else
    disp(['check on-surface nodes'])
    keyboard % check vioreanu-rokhlin nodes
  end
end
pt = s.thi - s.tlo;                  % panel size in parameter
t = zeros(1,N); s.w = t;

for i=1:np
    ii = (i-1)*p+(1:p); % indices of this panel
    t(ii) = s.tlo(i) + (1+x')/2*pt(i); s.w(ii) = w*pt(i)/2; % nodes weights this panel
end

if isfield(s,'Z')
  s.x = s.Z(t); % quadr nodes
else
  if x_uvs_qntype == 'V'
    if isfield(s,'x_mat')
      s.x = (s.x_mat*coefs_xyz)';
    else
      s.x = (koorn_coefs2vals_vec(order-1,order*(order+1)/2,[real(striZ(t));imag(striZ(t))],polsKN)*coefs_xyz)';
    end
  elseif x_uvs_qntype == 'T'
    if isfield(s,'x_mat')
      s.x = (s.x_mat*x_uvs')';
    else
      s.x = interpval([real(rboxZ(t));imag(rboxZ(t))],x_uvs,xtmp);
    end
  end
end

s.xp = zeros(size(s.x));
s.xpp = zeros(size(s.x));

if isfield(s,'Zp')
    s.xp = s.Zp(t);  % 1st dir of curve
else
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        s.xp(1,ii) = (D*s.x(1,ii)'*2/pt(i))'; s.xp(2,ii) = (D*s.x(2,ii)'*2/pt(i))'; s.xp(3,ii) = (D*s.x(3,ii)'*2/pt(i))';
    end
end
if isfield(s,'Zpp') 
    s.xpp = s.Zpp(t); % 2nd dir of curve
else
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        s.xpp(1,ii) = (D*s.xp(1,ii)'*2/pt(i))'; s.xpp(2,ii) = (D*s.xp(2,ii)'*2/pt(i))'; s.xpp(3,ii) = (D*s.xp(3,ii)'*2/pt(i))';
    end
end

s.sp = sqrt(s.xp(1,:).^2+s.xp(2,:).^2+s.xp(3,:).^2); 
s.tang = s.xp./(ones(3,1)*s.sp); 
s.w = s.w.*s.sp; % speed weights
s.t = t;
end

function test_quadr_3dline
close all
theta = 3;
s.Z = @(t) [(2*pi-t).*cos(theta*t)/(2*pi); (2*pi-t).*sin(theta*t)/(2*pi);t];
tt = linspace(0,2*pi,100);
temp = s.Z(tt);
figure(1), plot3(temp(1,:),temp(2,:),temp(3,:))

N = 200; qntype = 'C';
[s, N, np] = quadr_3dline(s, N, qntype);
figure(1), hold on, plot3(s.x(1,:),s.x(2,:),s.x(3,:),'.r');
quiver3(s.x(1,:),s.x(2,:),s.x(3,:),s.tang(1,:),s.tang(2,:),s.tang(3,:))

s.Zp = @(t) [(-(2*pi-t).*sin(theta*t)*theta-cos(theta*t))/(2*pi);
             ((2*pi-t).*cos(theta*t)*theta-sin(theta*t))/(2*pi);
             ones(size(t))];
test = s.Zp(s.t); 
diff = test-s.xp;
max(abs(diff(:)))

% arc_ds = @(t) sqrt(s.Zp{1}(t).^2+s.Zp{2}(t).^2+s.Zp{3}(t).^2);
% arclength = integral(arc_ds,0,2*pi);
% diff2 = arclength - sum(s.ws)

a = 2*pi; b = ((2*pi)^2+1)/theta^2;
true_arclength = @(t) abs(theta)/(4*pi)*((t-a)*sqrt((t-a)^2+b)-b*log(sqrt((t-a)^2+b)+a-t));
diff3 = sum(s.ws) - (true_arclength(2*pi)-true_arclength(0))

keyboard
end
