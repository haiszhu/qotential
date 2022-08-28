function y = transvector(nz,nx,y0)
% nz, normal at xc
ny = cross(nz,nx);  % nx, ny, nz, new coord x,y,z direction, also some of them are picked to be normal directions
y = [nx'*y0;ny'*y0;nz'*y0]; % projection on new coord
% y = y0 - ones(size(y0));
end