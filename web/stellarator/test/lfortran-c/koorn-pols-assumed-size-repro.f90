module koorn_pols_assumed_size_repro_mod
  implicit none
contains
  subroutine koorn_pols_assumed_size_repro(nmax, pols)
    integer(8), intent(in) :: nmax
    real(8), intent(out) :: pols((nmax+1_8)*(nmax+2_8)/2_8)
    real(8) :: legpols(0:100), jacpols(0:100,0:100)
    integer(8) :: n, k, iii

    legpols = 1.0_8
    jacpols = 1.0_8
    iii = 0_8
    do n = 0_8, nmax
      do k = 0_8, n
        iii = iii + 1_8
        pols(iii) = legpols(k) * jacpols(n-k,k)
      end do
    end do
  end subroutine koorn_pols_assumed_size_repro
end module koorn_pols_assumed_size_repro_mod
