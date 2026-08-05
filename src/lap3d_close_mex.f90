
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

subroutine qol_rrq_mex(m, tx, iside, order, h_dim, nquad, n, sx, sw, snx, &
                       rts, rps, orderff, distff, isimd, ichart, sxbd_chart, rv_chart, &
                       As, Ad, Omega, IalphaAsvestas, timeinfo)
  use lap3d_close_mod, only: rrq_r64
  implicit none
  integer(8), intent(in)    :: m, iside, order, h_dim, nquad, n
  integer(8), intent(in)    :: orderff, isimd
  real(8),    intent(in)    :: tx(3,m), sx(3,n), sw(n), snx(3,n)
  real(8),    intent(in)    :: rts(3,n), rps(3,n), distff
  real(8),    intent(inout) :: As(m,n), Ad(m,n)
  real(8),    intent(inout) :: Omega(4*h_dim,m), IalphaAsvestas(m)
  real(8),    intent(inout) :: timeinfo(20)
  integer(8), intent(in)    :: ichart
  real(8),    intent(in)    :: sxbd_chart(3,3*nquad), rv_chart(3,3)

  call rrq_r64(m, tx, n, sx, snx, sw, rts, rps, order, nquad, orderff, &
               distff, mod(iside,10_8) == 0_8, isimd, &
               ichart, sxbd_chart, rv_chart, &
               As, Ad, Omega, IalphaAsvestas, timeinfo)
end subroutine qol_rrq_mex
