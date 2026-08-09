module stellarator_wasm_api
  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int64_t
  use lap3d_close_mod, only: rrq_r64
  implicit none
  private
  public :: stellarator_rrq_r64

contains

  ! Narrow descriptor probe retained for LFortran backend regression tests.
  ! The full browser ABI is implemented by native/wasm_api_adapter.c around
  ! the generated stellarator_grf_core_mod C entry points.
  subroutine stellarator_rrq_r64(m, tx, n, sx, snx, sw, rts, rps, &
                                 order, nquad, orderff, distff, exterior, &
                                 isimd, ichart, sxbd_chart, rv_chart, &
                                 tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, &
                                 as, ad, omega, ialpha_asvestas, timeinfo) &
      bind(C, name="stellarator_rrq_r64")
    integer(c_int64_t), value, intent(in) :: m, n, order, nquad, orderff
    integer(c_int64_t), value, intent(in) :: isimd, ichart
    real(c_double), intent(in) :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(c_double), intent(in) :: rts(3,n), rps(3,n), distff
    logical(c_bool), value, intent(in) :: exterior
    real(c_double), intent(in) :: sxbd_chart(3,3*nquad), rv_chart(3,3)
    real(c_double), intent(in) :: tgl(nquad), wgl(nquad)
    real(c_double), intent(in) :: Dgl(nquad,nquad), w_bclag(nquad)
    real(c_double), intent(in) :: Legmat(nquad,nquad)
    real(c_double), intent(in) :: umatr(order*(order+1)/2,order*(order+1)/2)
    real(c_double), intent(in) :: vmatr(order*(order+1)/2,order*(order+1)/2)
    real(c_double), intent(out) :: as(m,n), ad(m,n)
    real(c_double), intent(out) :: omega(4*(order*(order+1)/2),m)
    real(c_double), intent(out) :: ialpha_asvestas(m)
    real(c_double), intent(inout) :: timeinfo(20)
    logical :: exterior_fortran

    exterior_fortran = exterior
    call rrq_r64(m, tx, n, sx, snx, sw, rts, rps, order, nquad, &
                 orderff, distff, exterior_fortran, isimd, ichart, sxbd_chart, &
                 rv_chart, tgl, wgl, Dgl, w_bclag, Legmat, umatr, vmatr, &
                 as, ad, omega, ialpha_asvestas, timeinfo)
  end subroutine stellarator_rrq_r64

end module stellarator_wasm_api
