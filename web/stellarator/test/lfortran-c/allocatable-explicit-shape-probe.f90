module allocatable_explicit_shape_probe_mod
  implicit none
contains
  subroutine caller(n, total)
    integer, intent(in) :: n
    real(8), intent(out) :: total
    real(8), allocatable :: values(:)
    allocate(values(n))
    call fill(n, values)
    total = sum(values)
  end subroutine caller

  subroutine fill(n, values)
    integer, intent(in) :: n
    real(8), intent(out) :: values(n)
    values = 1.0d0
  end subroutine fill
end module allocatable_explicit_shape_probe_mod
