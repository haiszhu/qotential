function [coeffs,powx,powy,powz,lens,ijidx]=hijcoeffsall_mex(order)

order2 = order*order;
mex_id_ = 'hijcoeffsall(i int[x], o double[xx], o int[xx], o int[xx], o int[xx], o int[x], o int[xx])';
[coeffs, powx, powy, powz, lens, ijidx] = specialquad(mex_id_, order, 1, order2, order2, order2, order2, order2, order2, order2, order2, order2, 2, order2);

end

