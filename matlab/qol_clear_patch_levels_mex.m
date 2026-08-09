function [ier] = qol_clear_patch_levels_mex(lv, ier)
lv = double(lv);
if nargin < 2 || isempty(ier), ier = 0; end
mex_id_ = 'qol_clear_patch_levels_mex(c i uint64_t[x], c io int64_t[x])';
[ier] = qotential_mex(mex_id_, lv, ier, 1, 1);
end
