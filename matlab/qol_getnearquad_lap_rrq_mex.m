function [sparsei, sparsej, sparsevs, sparsevd, timeinfo] = qol_getnearquad_lap_rrq_mex(npan, nterms, hdim, nquad, nsrc, sx, snx, sw, rts, rps, ntcx, tcx, tcxrow, tcxi, orderff, distff, iside, isimd, nnz, sparsei, sparsej, sparsevs, sparsevd, timeinfo)
npan = double(npan);
nterms = double(nterms);
hdim = double(hdim);
nquad = double(nquad);
nsrc = double(nsrc);
ntcx = double(ntcx);
orderff = double(orderff);
distff = double(distff);
iside = double(iside);
isimd = double(isimd);
nnz = double(nnz);
npanp1 = npan+1;
sx = double(reshape(sx,3,nsrc));
snx = double(reshape(snx,3,nsrc));
sw = double(reshape(sw,nsrc,1));
rts = double(reshape(rts,3,nsrc));
rps = double(reshape(rps,3,nsrc));
tcx = double(reshape(tcx,3,ntcx));
tcxrow = double(reshape(tcxrow,ntcx,1));
tcxi = double(reshape(tcxi,npanp1,1));
if hdim ~= nterms*(nterms+1)/2 || npan*hdim ~= nsrc
    error('qol_getnearquad_lap_rrq_mex:layout', ...
        'The panel order, hdim, npan, and nsrc are inconsistent.');
end
if any(tcxi ~= round(tcxi)) || tcxi(1) ~= 1 || ...
        any(diff(tcxi) < 0) || tcxi(end) ~= ntcx+1
    error('qol_getnearquad_lap_rrq_mex:layout', ...
        'tcxi must be monotone, one-based, and cover tcx exactly.');
end
required_nnz = (tcxi(end)-1)*hdim;
if required_nnz > nnz
    error('qol_getnearquad_lap_rrq_mex:sparseSize', ...
        'nnz is too small for the requested panel-target interactions.');
end
if nargin < 20 || isempty(sparsei), sparsei = zeros(nnz,1); end
if nargin < 21 || isempty(sparsej), sparsej = zeros(nnz,1); end
if nargin < 22 || isempty(sparsevs), sparsevs = zeros(nnz,1); end
if nargin < 23 || isempty(sparsevd), sparsevd = zeros(nnz,1); end
if nargin < 24 || isempty(timeinfo), timeinfo = zeros(20,1); end
mex_id_ = 'qol_getnearquad_lap_rrq_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[xx], c i double[xx], c i double[x], c i double[xx], c i double[xx], c i int64_t[x], c i double[xx], c i int64_t[x], c i int64_t[x], c i int64_t[x], c i double[x], c i int64_t[x], c i int64_t[x], c i int64_t[x], c io int64_t[x], c io int64_t[x], c io double[x], c io double[x], c io double[x])';
[sparsei, sparsej, sparsevs, sparsevd, timeinfo] = qotential_mex(mex_id_, npan, nterms, hdim, nquad, nsrc, sx, snx, sw, rts, rps, ntcx, tcx, tcxrow, tcxi, orderff, distff, iside, isimd, nnz, sparsei, sparsej, sparsevs, sparsevd, timeinfo, 1, 1, 1, 1, 1, 3, nsrc, 3, nsrc, nsrc, 3, nsrc, 3, nsrc, 1, 3, ntcx, ntcx, npanp1, 1, 1, 1, 1, 1, nnz, nnz, nnz, nnz, 20);
end

