function [Ac] = qol_lap3ddlp_closepanel_mex(t, s, order, ref, if_adapt)
t_x        = double(reshape(t.x, 3, []));
s_x        = double(reshape(s.x, 3, []));
m_tgt      = size(t_x, 2);
npat       = size(s_x, 2);
order      = double(order);
ref        = double(ref);
if_adapt_d = double(logical(if_adapt));
Ac         = zeros(m_tgt, npat);
mex_id_ = 'qol_lap3ddlp_closepanel_mex(c i int64_t[x], c i double[xx], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i double[x], c io double[xx])';
[Ac] = qotential_mex(mex_id_, m_tgt, t_x, npat, s_x, order, ref, if_adapt_d, Ac, 1, 3, m_tgt, 1, 3, npat, 1, 1, 1, m_tgt, npat);
end
