function [xx,yy] = meshgrid_mex(x,y)

n = numel(x);
m = numel(y);

mex_id_ = 'meshgrid(i int[x], i int[x], i int[x], i int[x], o int[xx], o int[xx])';
[xx, yy] = specialquad(mex_id_, n, x, m, y, 1, n, 1, m, m, n, m, n);

end

