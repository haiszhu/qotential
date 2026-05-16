function [xx,yy] = meshgrid_mex(x,y)

n = numel(x);
m = numel(y);

mex_id_ = 'meshgrid(c i int[x], c i int[x], c i int[x], c i int[x], c o int[xx], c o int[xx])';
[xx, yy] = specialquad(mex_id_, n, x, m, y, 1, n, 1, m, m, n, m, n);

end

