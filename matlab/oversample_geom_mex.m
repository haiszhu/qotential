function srcover = oversample_geom_mex(npatches,norders,ixyzs,iptype,npts,srccoefs,srcvals,novers,ixyzso,npts_over)

npatchesp1 = npatches+1;
srcover = zeros(12,npts_over);

mex_id_ = 'oversample_geom(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[xx])';
[srcover] = fmm3dbie_helpers_mex(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, srcvals, novers, ixyzso, npts_over, srcover, 1, npatches, npatchesp1, npatches, 1, 9, npts, 12, npts, npatches, npatchesp1, 1, 12, npts_over);

end

% ---------------------------------------------------------------------
