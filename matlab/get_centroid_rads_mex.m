function [cms,rads] = get_centroid_rads_mex(npatches,norders,ixyzs,iptype,npts,srccoefs)

npatchesp1 = npatches+1;
cms = zeros(3,npatches);
rads = zeros(npatches,1);

mex_id_ = 'get_centroid_rads(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c io double[xx], c io double[x])';
[cms, rads] = fmm3dbie_helpers_mex(mex_id_, npatches, norders, ixyzs, iptype, npts, srccoefs, cms, rads, 1, npatches, npatchesp1, npatches, 1, 9, npts, 3, npatches, npatches);

end
