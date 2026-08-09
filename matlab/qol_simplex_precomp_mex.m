function [tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr] = qol_simplex_precomp_mex(nquad, korder, kpols, tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr)
nquad = double(nquad);
korder = double(korder);
kpols = double(kpols);
if nargin < 4 || isempty(tgl), tgl = zeros(nquad, 1); end
if nargin < 5 || isempty(wgl), wgl = zeros(nquad, 1); end
if nargin < 6 || isempty(Dgl), Dgl = zeros(nquad, nquad); end
if nargin < 7 || isempty(w_bclag), w_bclag = zeros(nquad, 1); end
if nargin < 8 || isempty(Legmat), Legmat = zeros(nquad, nquad); end
if nargin < 9 || isempty(umatr), umatr = zeros(kpols, kpols); end
if nargin < 10 || isempty(vmatr), vmatr = zeros(kpols, kpols); end
mex_id_ = 'qol_simplex_precomp_mex(c i int64_t[x], c i int64_t[x], c i int64_t[x], c io double[x], c io double[x], c io double[xx], c io double[x], c io double[xx], c io double[xx], c io double[xx])';
[tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr] = qotential_mex(mex_id_, nquad, korder, kpols, tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, 1, 1, 1, nquad, nquad, nquad, nquad, nquad, nquad, nquad, kpols, kpols, kpols, kpols);
end
