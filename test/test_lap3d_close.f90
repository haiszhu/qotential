! qotential/test/test_lap3d_close.f90
!
! Smoke test for Lap3dDLP_closepanel_r64.  Builds a small torus-patch
! source surface (same shape demo1.m uses, order=8) and a handful of
! close target points, calls Lap3dDLP_closepanel_r64 on the adaptive
! path, and reports basic sanity of the returned Ac matrix:
!
!   - shape is (m_tgt, npat)
!   - no NaNs / Infs
!   - row magnitudes are in a sensible range
!
! Not a numerical-accuracy test.  For that, compare against the MATLAB
! Lap3dDLP_closepanel_demo via the mex wrappers (separate task).

program test_lap3d_close
  use lap3d_close_mod, only: Lap3dDLP_closepanel_r64
  implicit none
  integer, parameter :: r64 = 8

  integer(8) :: order, ref, npat, m_tgt
  integer(8) :: i, j, k, ip
  real(r64), allocatable :: s_x(:,:), t_x(:,:), Ac(:,:)
  real(r64) :: tgl_v(8), wgl_v(8), x1, x2
  real(r64) :: a_torus, b_torus, dt, dp, tc, pc
  real(r64) :: rmean(3), max_abs, min_abs
  logical   :: any_nonfinite

  order = 8_8
  ref   = 1_8
  npat  = order * order
  m_tgt = 5_8

  ! ---- build a small torus-patch (a=1, b=0.5, panel-sized centered) ----
  ! GL nodes for order=8 (hard-coded so we don't need to call gauss here)
  tgl_v(1) = -0.9602898564975363_r64; wgl_v(1) = 0.1012285362903763_r64
  tgl_v(2) = -0.7966664774136267_r64; wgl_v(2) = 0.2223810344533745_r64
  tgl_v(3) = -0.5255324099163290_r64; wgl_v(3) = 0.3137066458778873_r64
  tgl_v(4) = -0.1834346424956498_r64; wgl_v(4) = 0.3626837833783620_r64
  tgl_v(5) =  0.1834346424956498_r64; wgl_v(5) = 0.3626837833783620_r64
  tgl_v(6) =  0.5255324099163290_r64; wgl_v(6) = 0.3137066458778873_r64
  tgl_v(7) =  0.7966664774136267_r64; wgl_v(7) = 0.2223810344533745_r64
  tgl_v(8) =  0.9602898564975363_r64; wgl_v(8) = 0.1012285362903763_r64

  a_torus = 1.0_r64
  b_torus = 0.5_r64
  dt = 0.5_r64 / real(ref, r64)     ! tpansiz1
  dp = 0.25_r64 / real(ref, r64)    ! tpansiz2
  tc = 0.0_r64
  pc = 0.0_r64

  allocate(s_x(3, npat))
  ip = 0_8
  do i = 1, order
    do j = 1, order
      ip = ip + 1_8
      x1 = tc + dt * tgl_v(j)        ! note: MATLAB orig+tpansiz1*x1(:)' indexing
      x2 = pc + dp * tgl_v(i)
      s_x(1, ip) = b_torus*sin(x1)
      s_x(2, ip) = (a_torus + b_torus*cos(x1)) * sin(x2)
      s_x(3, ip) = (a_torus + b_torus*cos(x1)) * cos(x2)
    end do
  end do
  ! Center the patch on origin (matches demo1's "r = r - mean(r,2)").
  rmean = 0.0_r64
  do ip = 1, npat
    rmean = rmean + s_x(:, ip)
  end do
  rmean = rmean / real(npat, r64)
  do ip = 1, npat
    s_x(:, ip) = s_x(:, ip) - rmean
  end do

  ! ---- a few targets near the patch centroid --------------------------
  allocate(t_x(3, m_tgt))
  do k = 1, m_tgt
    t_x(1, k) = 0.05_r64 * real(k, r64)
    t_x(2, k) = 0.03_r64 * real(k, r64) - 0.05_r64
    t_x(3, k) = 0.02_r64 * real(k, r64) + 0.01_r64
  end do

  ! ---- call Lap3dDLP_closepanel_r64 -----------------------------------
  allocate(Ac(m_tgt, npat))
  Ac = 0.0_r64
  call Lap3dDLP_closepanel_r64(m_tgt, t_x, npat, s_x, order, ref, .true., Ac)

  ! ---- sanity checks --------------------------------------------------
  any_nonfinite = .false.
  max_abs       = 0.0_r64
  min_abs       = huge(1.0_r64)
  do k = 1, m_tgt
    do ip = 1, npat
      if (.not. (abs(Ac(k, ip)) < huge(1.0_r64))) then
        any_nonfinite = .true.
      end if
      if (abs(Ac(k, ip)) > max_abs) max_abs = abs(Ac(k, ip))
      if (abs(Ac(k, ip)) < min_abs) min_abs = abs(Ac(k, ip))
    end do
  end do

  write(*,'(A)')       '=== test_lap3d_close ==='
  write(*,'(A,I0,A,I0)') '  order = ', order, '   ref = ', ref
  write(*,'(A,I0,A,I0)') '  m_tgt = ', m_tgt, '   npat = ', npat
  write(*,'(A,ES12.4)')  '  max|Ac| = ', max_abs
  write(*,'(A,ES12.4)')  '  min|Ac| = ', min_abs
  if (any_nonfinite) then
    write(*,'(A)') '  FAIL: Ac contains non-finite entries'
    error stop
  else
    write(*,'(A)') '  PASS: Ac shape correct, all entries finite.'
  end if

end program test_lap3d_close
