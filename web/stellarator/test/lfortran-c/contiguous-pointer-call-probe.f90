module contiguous_pointer_call_probe_mod
  implicit none
contains
  subroutine fill_values(n, values)
    integer(8), intent(in) :: n
    real(8), intent(out) :: values(n)
    values = 1.0_8
  end subroutine fill_values

  subroutine pointer_caller(n, total)
    integer(8), intent(in) :: n
    real(8), intent(out) :: total
    real(8), pointer, contiguous :: values(:)

    allocate(values(n))
    call fill_values(n, values)
    total = sum(values)
    deallocate(values)
  end subroutine pointer_caller
end module contiguous_pointer_call_probe_mod
