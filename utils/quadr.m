function [s, N, np] = quadr(s, N, qtype, qntype, x, w, D)  % set up quadrature on a closed segment
% QUADR - set up quadrature (either global or panel-based) on a segment struct
%
% [s N] = quadr(s, N, qtype, qntype) adds quadrature input to segment struct s.
% Inputs:
%  s - segment struct containing parametrization
%  s.Z - complex function, parameterization of the segment on [0,2\pi]
%  s.T (optional) - reparameterization function mapping [0,2\pi] to itself.
%  N - requested number of nodes
%  qtype - quadrature type: 'g' global periodic trapezoid rule
%                           'p' panel-wise quadr w/ s.p nodes per pan
%  qntype - type of panel quadr: 'G' Gause-Legendre
%                                'C' Cheyshev (Fejer's 1st rule)
%
% Outputs: s - amended segment struct
%          N - actual number of nodes (is within s.p of requested N)
%          np - number of panels
% Notes: 1) Note sign change in normal sense vs periodicdirpipe.m
% 2) Non-adaptive for now.  Barnett 6/10/13
% - Added reparameterized panel. Bowei Wu 3/22/18

if nargin < 3
    qtype = 'g'; 
elseif nargin < 4 && qtype=='p'
    qntype = 'G'; 
else
    p = s.p; if qntype=='G', [x, w, D] = gauss(p); else, [x, w, D] = cheby(p); end  
end

if qtype=='g' % global (ptr) quadr
    if isfield(s,'Z') && ~isempty(N) % quadr nodes
        t = (1:N)'/N*2*pi;
        s.x = s.Z(t);
        s.tlo = 0; s.thi = 2*pi; s.p = N; s.w = 2*pi/N*ones(N,1); np=1; % 1 big panel
        s.xlo = s.Z(s.tlo); s.xhi = s.Z(s.thi);
    elseif isfield(s,'x')
        s.x = s.x(:);
    else
        error('Need to provide at least s.Z and N, or s.x. Neither found!');
    end
elseif qtype=='p'
    if ~isfield(s,'Z'), error('Need to provide s.Z to build panels!'); end 
    if ~isfield(s,'p'), s.p=16; end       % default panel order
    p = s.p;
    np = 1;          % default, to prevent failure of line s.np=np below 4/4/19
    if isfield(s,'tpan')   % adaptive panel, See adaptive_panel.m
        s.tlo = s.tpan(1:end-1);
        s.thi = s.tpan(2:end);
        np = numel(s.tlo);
        N = np*p;
    elseif ~isfield(s,'tlo') || ~isfield(s,'thi') % if panels are not given
      if isempty(N), error('Need to provide N (approx num of pts) to build panels'), end
        np = ceil(N/p); N = p*np;      % np = # panels
        s.tlo = (0:np-1)'/np*2*pi; % panel start params
        s.thi = (1:np)'/np*2*pi; % panel end params
        if isfield(s,'T')   % adaptive panel via reparameterization, See reparam_mixed.m
            s.tlo = s.T(s.tlo);
            s.thi = s.T(s.thi);
        end
    end
    s.np = np;
    s.xlo = s.Z(s.tlo); % panel start locs
    s.xhi = s.Z(s.thi);  % panel end locs
    pt = s.thi - s.tlo;                  % panel size in parameter
    t = zeros(N,1); s.w = t;
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        t(ii) = s.tlo(i) + (1+x)/2*pt(i); s.w(ii) = w*pt(i)/2; % nodes weights this panel
    end
    s.x = s.Z(t); % quadr nodes
end

if N~=length(s.x), N = length(s.x); warning('N has changed!'); end 
s.xp = zeros(length(s.x),1);
s.xpp = zeros(length(s.x),1);

if isfield(s,'Zp'), s.xp = s.Zp(t); % 1st dir of curve
else
    if qtype == 'p'
        for i=1:np
            ii = (i-1)*p+(1:p); % indices of this panel
            s.xp(ii) = D*s.x(ii)*2/pt(i);
        end
    else
        s.xp = perispecdiff(s.x); % fourier spectral diff
    end
end
if isfield(s,'Zpp'), s.xpp = s.Zpp(t); % 2nd dir of curve
else
    if qtype == 'p'
    for i=1:np
        ii = (i-1)*p+(1:p); % indices of this panel
        s.xpp(ii) = D*s.xp(ii)*2/pt(i);
    end
    else
        s.xpp = perispecdiff(s.xp); % fourier spectral diff
    end
end

s.sp = abs(s.xp); s.tang = s.xp./s.sp; s.nx = -1i*s.tang; % speed, tangent, normal
s.cur = -real(conj(s.xpp).*s.nx)./s.sp.^2; % curvature
s.ws = s.w.*s.sp; % speed weights
s.t = t; s.wxp = s.w.*s.xp; % complex speed weights (Helsing's wzp)
end

function g = perispecdiff(f)
% PERISPECDIFF - use FFT to take periodic spectral differentiation of vector
%
% g = PERISPECDIFF(f) returns g the derivative of the spectral interpolant
%  of f, which is assumed to be the values of a smooth 2pi-periodic function
%  at the N gridpoints 2.pi.j/N, for j=1,..,N (or any translation of such
%  points).
%
% Barnett 2/18/14
N = numel(f);
if mod(N,2)==0   % even
  g = ifft(fft(f(:)).*[0 1i*(1:N/2-1) 0 1i*(-N/2+1:-1)].');
else
  g = ifft(fft(f(:)).*[0 1i*(1:(N-1)/2) 1i*((1-N)/2:-1)].');
end
g = reshape(g,size(f));
end
