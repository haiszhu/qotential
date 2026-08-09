function [ier] = qol_clear_patch_levels_all_mex(ier)
if nargin < 1 || isempty(ier), ier = 0; end
mex_id_ = 'qol_clear_patch_levels_all_mex(c io int64_t[x])';
[ier] = qotential_mex(mex_id_, ier, 1);
end
