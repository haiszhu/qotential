function [norders,ixyzs,iptype,srcvals,srccoefs,wts] = get_wtorus_geom_mex(radii,scales,nosc,nu,nv,norder)

npatches = 2*nu*nv;
npols = (norder+1)*(norder+2)/2;
npts = npatches*npols;
npatchesp1 = npatches+1;

norders = zeros(npatches,1);
ixyzs = zeros(npatchesp1,1);
iptype = zeros(npatches,1);
srcvals = zeros(12,npts);
srccoefs = zeros(9,npts);
wts = zeros(npts,1);

mex_id_ = 'get_wtorus_geom(c i double[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io int64_t[x], c io double[xx], c io double[xx], c io double[x])';
[norders, ixyzs, iptype, srcvals, srccoefs, wts] = fmm3dbie_helpers_mex(mex_id_, radii, scales, nosc, nu, nv, npatches, norder, npts, norders, ixyzs, iptype, srcvals, srccoefs, wts, 3, 3, 1, 1, 1, 1, 1, 1, npatches, npatchesp1, npatches, 12, npts, 9, npts, npts);

end

% ---------------------------------------------------------------------
