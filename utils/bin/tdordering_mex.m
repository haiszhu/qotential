function  ijidx = tdordering_mex(orderx,ordery)

orderxy = orderx*ordery;
mex_id_ = 'tdordering(i int[x], i int[x], o int[xx])';
[ijidx] = specialquad(mex_id_, orderx, ordery, 1, 1, 2, orderxy);

end

