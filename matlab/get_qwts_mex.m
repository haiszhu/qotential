function wover = get_qwts_mex(npatches,novers,ixyzso,iptype,npts_over,srcover)

npatchesp1 = npatches+1;
wover = zeros(npts_over,1);

mex_id_ = 'get_qwts(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c io double[x])';
[wover] = fmm3dbie_helpers_mex(mex_id_, npatches, novers, ixyzso, iptype, npts_over, srcover, wover, 1, npatches, npatchesp1, npatches, 1, 12, npts_over, npts_over);

end

% ---------------------------------------------------------------------
