function y = transcoord(xc,nz,nx,origin,y0)
% nz, normal at xc
xc = xc+origin*nz;
ny = cross(nz,nx);  % nx, ny, nz, new coord x,y,z direction
y = [nx'*(y0-xc*ones(1,size(y0,2)));ny'*(y0-xc*ones(1,size(y0,2)));nz'*(y0-xc*ones(1,size(y0,2)))];
% y = y0 - ones(size(y0));