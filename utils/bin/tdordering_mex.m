function  ijidx = tdordering_mex(orderx,ordery)

orderxy = orderx*ordery;
mex_id_ = 'tdordering(c i int[x], c i int[x], c o int[xx])';
[ijidx] = specialquad(mex_id_, orderx, ordery, 1, 1, 2, orderxy);

end

