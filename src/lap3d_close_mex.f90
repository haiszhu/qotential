! qotential/src/lap3d_close_mex.f90
!
! Top-level mwrap wrapper(s) for the BIE-solver glue routines in
! lap3d_close_mod.  Symbol prefix: qol_ (QOtential, Lap3d-close).
!
! No module-mangling -- these are free-standing externals so mwrap
! can bind to them directly via the # FORTRAN line in qotential.mw.
! All arguments cross the boundary as integer(8) or double precision
! (matlab-fortran-skill conventions).

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
