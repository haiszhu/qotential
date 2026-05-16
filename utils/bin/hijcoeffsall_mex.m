function [coeffs,powx,powy,powz,lens,ijidx]=hijcoeffsall_mex(order)

order2 = order*order;
mex_id_ = 'hijcoeffsall(c i int[x], c o double[xx], c o int[xx], c o int[xx], c o int[xx], c o int[x], c o int[xx])';
[coeffs, powx, powy, powz, lens, ijidx] = specialquad(mex_id_, order, 1, order2, order2, order2, order2, order2, order2, order2, order2, order2, 2, order2);

end

