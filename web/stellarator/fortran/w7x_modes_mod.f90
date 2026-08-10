! Generated from test/w7x_grf.m by generate-w7x-modes.py.
! Do not edit this coefficient table by hand.
module w7x_modes_mod
  use quatapproximation_mod, only: r64
  implicit none
  private
  integer(8), parameter, public :: W7X_NMODE = 288_8
  integer(8), parameter, public :: W7X_NFP = 5_8
  public :: load_w7x_modes

contains

  subroutine load_w7x_modes(mn, rc, zs)
    real(r64), intent(out) :: mn(2*W7X_NMODE), rc(W7X_NMODE), zs(W7X_NMODE)
    include 'w7x-modes-dat.txt'
  end subroutine load_w7x_modes
end module w7x_modes_mod
