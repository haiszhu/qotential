function [lv, ier] = qol_build_patch_levels_mex(korder, npols, sx, rts, rps, order, dist, nlevel, lv, ier)
korder = double(korder);
npols = double(npols);
order = double(order);
nlevel = double(nlevel);
dist = double(dist);
sx = double(reshape(sx, 3, npols));
rts = double(reshape(rts, 3, npols));
rps = double(reshape(rps, 3, npols));
if nargin < 9 || isempty(lv), lv = 0; end
if nargin < 10 || isempty(ier), ier = 0; end
mex_id_ = 'qol_build_patch_levels_mex(c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[xx], c i int64_t[x], c i double[x], c i int64_t[x], c io uint64_t[x], c io int64_t[x])';
[lv, ier] = qotential_mex(mex_id_, korder, npols, sx, rts, rps, order, dist, nlevel, lv, ier, 1, 1, 3, npols, 3, npols, 3, npols, 1, 1, 1, 1, 1);
end
