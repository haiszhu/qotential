
subroutine qol_lap3ddlp_closepanel_mex(m_tgt, t_x, npat, s_x, order, ref, &
                                        if_adapt_d, Ac)
  use lap3d_close_mod, only: Lap3dDLP_closepanel_r64
  implicit none
  integer(8), intent(in)    :: m_tgt, npat, order, ref
  real(8),    intent(in)    :: t_x(3, m_tgt)
  real(8),    intent(in)    :: s_x(3, npat)
  real(8),    intent(in)    :: if_adapt_d           ! 0.0 = false, nonzero = true
  real(8),    intent(inout) :: Ac(m_tgt, npat)

  logical :: if_adapt
  if_adapt = (if_adapt_d /= 0.0d0)
  call Lap3dDLP_closepanel_r64(m_tgt, t_x, npat, s_x, order, ref, if_adapt, Ac)
end subroutine qol_lap3ddlp_closepanel_mex


subroutine qol_rrq_mex(m, tx, n, sx, snx, sw, rts, rps, order, nquad, &
                       orderff, distff, exterior, isimd, ichart, &
                       sxbd_chart, rv_chart, tgl, wgl, Dgl, w_bclag, &
                       Legmat, umatr, vmatr, As, Ad, Omega, &
                       IalphaAsvestas, timeinfo)
  use lap3d_close_mod, only: rrq_r64
  implicit none
  integer(8), intent(in)    :: m, n, order, nquad, orderff, isimd
  real(8),    intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
  real(8),    intent(in)    :: rts(3,n), rps(3,n), distff
  logical,    intent(in)    :: exterior
  real(8),    intent(out)   :: As(m,n), Ad(m,n)
  real(8),    intent(out)   :: Omega(4*(order*(order+1)/2),m)
  real(8),    intent(out)   :: IalphaAsvestas(m)
  real(8),    intent(inout) :: timeinfo(20)
  integer(8), intent(in)    :: ichart
  real(8),    intent(in)    :: sxbd_chart(3,3*nquad), rv_chart(3,3)
  real(8),    intent(in)    :: tgl(nquad), wgl(nquad)
  real(8),    intent(in)    :: Dgl(nquad,nquad), w_bclag(nquad)
  real(8),    intent(in)    :: Legmat(nquad,nquad)
  real(8),    intent(in)    :: umatr(order*(order+1)/2,order*(order+1)/2)
  real(8),    intent(in)    :: vmatr(order*(order+1)/2,order*(order+1)/2)

  call rrq_r64(m, tx, n, sx, snx, sw, rts, rps, order, nquad, orderff, &
               distff, exterior, isimd, ichart, sxbd_chart, rv_chart, &
               tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, &
               As, Ad, Omega, IalphaAsvestas, timeinfo)
end subroutine qol_rrq_mex


subroutine qol_simplex_precomp_mex(nquad, korder, kpols, tgl, wgl, Dgl, &
    w_bclag, Legmat, umatr, vmatr)
  use lap3d_close_mod, only: simplex_precomp_r64
  implicit none
  integer(8), intent(in) :: nquad, korder, kpols
  real(8), intent(out) :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8), intent(out) :: w_bclag(nquad), Legmat(nquad,nquad)
  real(8), intent(out) :: umatr(kpols,kpols), vmatr(kpols,kpols)

  call simplex_precomp_r64(nquad, korder, kpols, tgl, wgl, Dgl, &
      w_bclag, Legmat, umatr, vmatr)
end subroutine qol_simplex_precomp_mex


subroutine qol_getnearquad_lap_rrq_mex(npan, nterms, hdim, nquad, &
    nsrc, sx, snx, sw, rts, rps, ntcx, tcx, tcxrow, tcxi, &
    orderff, distff, iside, isimd, nnz, sparsei, sparsej, &
    sparsevs, sparsevd, timeinfo)
  use lap3d_close_mod, only: getnearquad_lap_rrq_r64
  implicit none
  integer(8), intent(in) :: npan, nterms, hdim, nquad, nsrc, ntcx
  integer(8), intent(in) :: orderff, iside, isimd, nnz
  real(8), intent(in) :: sx(3,nsrc), snx(3,nsrc), sw(nsrc)
  real(8), intent(in) :: rts(3,nsrc), rps(3,nsrc)
  real(8), intent(in) :: tcx(3,ntcx), distff
  integer(8), intent(in) :: tcxrow(ntcx), tcxi(npan+1)
  integer(8), intent(inout) :: sparsei(nnz), sparsej(nnz)
  real(8), intent(inout) :: sparsevs(nnz), sparsevd(nnz), timeinfo(20)

  call getnearquad_lap_rrq_r64(npan, nterms, hdim, nquad, &
      nsrc, sx, snx, sw, rts, rps, ntcx, tcx, tcxrow, tcxi, &
      orderff, distff, iside, isimd, nnz, sparsei, sparsej, &
      sparsevs, sparsevd, timeinfo)
end subroutine qol_getnearquad_lap_rrq_mex


subroutine qol_build_closepanel_precomp_mex(n, sx, snx, sw, r_vert, &
    nterms, h_dim, nbd, sbdnp, nquad, alpha, exterior, ichart, sxbd_chart, &
    tgl, wgl, Dgl, umatr, R, c, sxbd, swbd, stangbd, sspbd, &
    qhat, Fbd, Fxbd, Fybd, Fzbd, Mmatrix)
  use lap3d_close_mod, only: build_closepanel_precomp_r64
  implicit none
  integer(8), intent(in) :: n, nterms, h_dim, nbd, sbdnp, nquad
  logical, intent(in) :: exterior
  integer(8), intent(in) :: ichart
  real(8), intent(in) :: sx(3,n), snx(3,n), sw(n), r_vert(3,3), alpha
  real(8), intent(in) :: sxbd_chart(3,nbd)
  real(8), intent(in) :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
  real(8), intent(in) :: umatr(h_dim,h_dim)
  real(8), intent(out) :: R(3,3), c(3), sxbd(3,nbd), swbd(nbd)
  real(8), intent(out) :: stangbd(3,nbd), sspbd(nbd), qhat(3)
  complex(8), intent(out) :: Fbd(nbd,(nterms+1)**2)
  complex(8), intent(out) :: Fxbd(nbd,(nterms+1)**2)
  complex(8), intent(out) :: Fybd(nbd,(nterms+1)**2)
  complex(8), intent(out) :: Fzbd(nbd,(nterms+1)**2)
  real(8), intent(out) :: Mmatrix(4*h_dim,4*h_dim)

  call build_closepanel_precomp_r64(n, sx, snx, sw, r_vert, nterms, &
      h_dim, nbd, sbdnp, nquad, alpha, exterior, ichart, &
      sxbd_chart, tgl, wgl, Dgl, umatr, R, c, sxbd, swbd, &
      stangbd, sspbd, qhat, Fbd, Fxbd, Fybd, Fzbd, Mmatrix)
end subroutine qol_build_closepanel_precomp_mex


subroutine qol_build_tc_chialpha_mex(m, r0, nbd, sbdnp, nquad, ncoeff, &
    sxbd, stangbd, sspbd, w1, w3, Fx, Fy, Fz, &
    IalphaAsvestas, Ichi, Ialpha)
  use lap3d_close_mod, only: build_tc_chialpha_r64
  implicit none
  integer(8), intent(in) :: m, nbd, sbdnp, nquad, ncoeff
  real(8), intent(in) :: r0(3,m), sxbd(3,nbd), stangbd(3,nbd)
  real(8), intent(in) :: sspbd(nbd)
  real(8), intent(in) :: w1(nquad,sbdnp,m), w3(nquad,sbdnp,m)
  complex(8), intent(in) :: Fx(nbd,ncoeff), Fy(nbd,ncoeff), Fz(nbd,ncoeff)
  real(8), intent(in) :: IalphaAsvestas(m)
  complex(8), intent(out) :: Ichi(m,ncoeff,4), Ialpha(m,ncoeff,4)
  call build_tc_chialpha_r64(m, r0, nbd, sbdnp, nquad, ncoeff, sxbd, &
      stangbd, sspbd, w1, w3, Fx, Fy, Fz, IalphaAsvestas, Ichi, Ialpha)
end subroutine qol_build_tc_chialpha_mex
