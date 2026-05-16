function [gradf,ijidx] = evaltensorproductharmonicgrad_mex(r,order)

gradf = [];
nt = numel(r(1,:));
orderp1 = order + 1;
order2p1 = 2*order + 1;
order2 = order^2;
mex_id_ = 'evaltensorproductharmonicgrad(c i int[x], c i double[xx], c i int[x], c o double[xx], c o double[xx], c o double[xx], c o double[xx], c o int[xx])';
[fx, fy, fz, f, ijidx] = specialquad(mex_id_, nt, r, order, 1, 3, nt, 1, nt, order2, nt, order2, nt, order2, nt, order2, 2, order2);

gradf.F1 = fx; gradf.F2 = fy; gradf.F3 = fz;
gradf.F = f;

% keyboard

end

