function [As, Ad, idxs, ms, ier] = qol_lap3dsdlpmat_levels_mex(m, tx, npols, lv, isimd, As, Ad, idxs, ms, ier)
m = double(m);
npols = double(npols);
lv = double(lv);
isimd = double(isimd);
tx = double(reshape(tx, 3, m));
if nargin < 6 || isempty(As), As = zeros(m, npols); end
if nargin < 7 || isempty(Ad), Ad = zeros(m, npols); end
if nargin < 8 || isempty(idxs), idxs = zeros(m, 1); end
if nargin < 9 || isempty(ms), ms = 0; end
if nargin < 10 || isempty(ier), ier = 0; end
mex_id_ = 'qol_lap3dsdlpmat_levels_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i uint64_t[x], c i int64_t[x], c io double[xx], c io double[xx], c io int64_t[x], c io int64_t[x], c io int64_t[x])';
[As, Ad, idxs, ms, ier] = qotential_mex(mex_id_, m, tx, npols, lv, isimd, As, Ad, idxs, ms, ier, 1, 3, m, 1, 1, 1, m, npols, m, npols, m, 1, 1);
end
