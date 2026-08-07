program array_is_contiguous_c
    implicit none
    real(8), allocatable :: values(:), strided(:)

    allocate(values(4))
    call fill_values(values)
    if (any(values /= 3.0_8)) error stop "contiguous"

    allocate(strided(8))
    strided = -1.0_8
    call fill_values(strided(1:8:2))
    if (any(strided(1:8:2) /= 3.0_8)) error stop "strided output"
    if (any(strided(2:8:2) /= -1.0_8)) error stop "strided untouched"
    deallocate(values)
    deallocate(strided)
contains
    subroutine fill_values(output)
        real(8), intent(out) :: output(4)
        output = 3.0_8
    end subroutine fill_values
end program array_is_contiguous_c
