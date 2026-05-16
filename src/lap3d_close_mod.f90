! qotential/src/lap3d_close_mod.f90
!
! BIE-solver layer: Fortran orchestration of the close-panel evaluation
! pipeline.  Glues the QuatApproximation-legacy package (close-eval
! kernel-coefficient builders, harmonic basis, boundary geometry) with
! the LineQuaaadrature-legacy package (moments / adaptive line quadrature)
! into a single per-patch routine.
!
! Mirrors the MATLAB function `Lap3dDLP_closepanel_demo` (local helper
! in qotential/demo1.m); produces the same close-evaluation matrix Ac.
!
! Public:
!   Lap3dDLP_closepanel_r64
!   Lap3dDLP_closepanel_r128   (skeleton; see r128 gap notes inside)
!
module lap3d_close_mod

  ! kinds + small primitives
  use quatapproximation_mod, only: r64, r128, gauss_r64, gauss_r128
  ! harmonic basis values and gradients
  use harmonic_mod, only: evaltensorproductharmonicgrad_r64
  ! boundary line quadrature on the tensor reference square
  use tensor_geom_mod, only: line3quadr_3dline_T, line3quadr_3dline_T_r128
  ! q^{nm} kernel-coefficient builders
  use qkernel_mod, only: qak_qnm_i_r64, qak_qnm_i_r128,  &
                          QAK_LPTYPE_D
  ! omega assembly (per-node and target-dependent)
  use omega_mod, only: qao_omeganm_i_r64, qao_omeganm_i_r128,         &
                       qao_omegaall_r64, qao_omegaall_r128
  ! in-house small LU solve (no LAPACK dep)
  use koorn_geom_mod, only: lu_solve_r128
  ! moments (from LineQuaaadrature-legacy)
  use solidangle_mod, only: eval_moments_funvals_r64

  implicit none
  private
  public :: Lap3dDLP_closepanel_r64
  public :: Lap3dDLP_closepanel_r128

  ! Geometry constant: the unit square boundary is traversed in 8
  ! equal panels (matches qotential demo1's linspace(0, 2*pi, 9)).
  integer(8), parameter :: SBDNP = 8_8

contains

  ! ================================================================
  ! Lap3dDLP_closepanel_r64
  !
  ! Build the close-evaluation Laplace DLP matrix Ac for one tensor-
  ! product patch and a set of close targets.  Single-precision r64.
  !
  ! Inputs:
  !   m_tgt              number of targets
  !   t_x(3, m_tgt)      target points
  !   npat               source-patch nodes count (= order^2)
  !   s_x(3, npat)       source-patch values at GL tensor grid
  !   order              patch tensor order  (sqrt(npat))
  !   ref                panel refinement factor (nquad = ref * order)
  !   if_adapt           .true. -> adaptive line-quad via
  !                      eval_moments_funvals_r64; .false. -> plain
  !                      smooth quadrature path (NOT YET PORTED here,
  !                      will error_stop -- use the if_adapt=.true.
  !                      path or call qotential's MATLAB
  !                      momentsallplain in the meantime).
  ! Output:
  !   Ac(m_tgt, npat)    close-evaluation matrix
  ! ================================================================
  subroutine Lap3dDLP_closepanel_r64(m_tgt, t_x, npat, s_x, order, ref, &
                                     if_adapt, Ac)
    integer(8), intent(in)  :: m_tgt, npat, order, ref
    real(r64),  intent(in)  :: t_x(3, m_tgt)
    real(r64),  intent(in)  :: s_x(3, npat)
    logical,    intent(in)  :: if_adapt
    real(r64),  intent(out) :: Ac(m_tgt, npat)

    integer(8) :: nquad, nbd, h_dim, morder, ncol, idx, k, kk
    real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:), tpan(:)
    real(r64), allocatable :: sxbd(:,:), swbd(:), stangbd(:,:), sspbd(:)
    real(r64), allocatable :: F(:,:), Fx(:,:), Fy(:,:), Fz(:,:)
    real(r64), allocatable :: F0(:,:), F1(:,:), F2(:,:), F3(:,:)
    integer(8), allocatable :: ijidx(:,:)
    real(r64), allocatable :: M_all(:,:,:)
    real(r64), allocatable :: q_i(:,:,:), q_j(:,:,:), q_k(:,:,:)
    real(r64), allocatable :: onm0(:,:,:), onm1(:,:,:), onm2(:,:,:), onm3(:,:,:)
    real(r64), allocatable :: dr(:,:)
    real(r64), allocatable :: Omega_all(:,:)
    real(r64), allocatable :: Mmat(:,:), rhs(:,:)
    real(r64), parameter :: PI = 4.0_r64*atan(1.0_r64)

    ! -- 1. quadrature setup ------------------------------------------
    nquad = ref * order
    nbd   = SBDNP * nquad
    h_dim = order * order
    morder = 2_8*order + 2_8           ! ncol for the M-block
    ncol  = morder

    allocate(tgl(nquad), wgl(nquad), Dgl(nquad, nquad))
    call gauss_r64(nquad, tgl, wgl, Dgl)

    allocate(tpan(SBDNP+1))
    do idx = 1, SBDNP+1
      tpan(idx) = 2.0_r64*PI * real(idx-1, r64) / real(SBDNP, r64)
    end do

    allocate(sxbd(3, nbd), swbd(nbd), stangbd(3, nbd), sspbd(nbd))
    call line3quadr_3dline_T(s_x, order, nquad, tgl, wgl, Dgl, &
                             SBDNP, tpan, nbd, sxbd, swbd, stangbd, sspbd)

    ! -- 2. harmonic gradient at source nodes -------------------------
    allocate(F(npat, h_dim), Fx(npat, h_dim), Fy(npat, h_dim), Fz(npat, h_dim))
    allocate(ijidx(2, h_dim))
    F = 0.0_r64; Fx = 0.0_r64; Fy = 0.0_r64; Fz = 0.0_r64; ijidx = 0_8
    call evaltensorproductharmonicgrad_r64(npat, s_x, order, Fx, Fy, Fz, F, ijidx)

    allocate(F0(h_dim, h_dim), F1(h_dim, h_dim), F2(h_dim, h_dim), F3(h_dim, h_dim))
    F0 = 0.0_r64
    F1 = Fx ; F2 = Fy ; F3 = Fz   ! 1st h_dim source-nodes (assumes npat=h_dim)

    ! Mmat = [ F0 -F1 -F2 -F3 ;
    !          F1  F0 -F3  F2 ;
    !          F2  F3  F0 -F1 ;
    !          F3 -F2  F1  F0 ]   (4*h_dim x 4*h_dim)
    allocate(Mmat(4*h_dim, 4*h_dim))
    Mmat = 0.0_r64
    Mmat(           1:  h_dim,           1:  h_dim) =  F0
    Mmat(           1:  h_dim,   h_dim+1:2*h_dim) = -F1
    Mmat(           1:  h_dim, 2*h_dim+1:3*h_dim) = -F2
    Mmat(           1:  h_dim, 3*h_dim+1:4*h_dim) = -F3
    Mmat(   h_dim+1:2*h_dim,           1:  h_dim) =  F1
    Mmat(   h_dim+1:2*h_dim,   h_dim+1:2*h_dim) =  F0
    Mmat(   h_dim+1:2*h_dim, 2*h_dim+1:3*h_dim) = -F3
    Mmat(   h_dim+1:2*h_dim, 3*h_dim+1:4*h_dim) =  F2
    Mmat( 2*h_dim+1:3*h_dim,           1:  h_dim) =  F2
    Mmat( 2*h_dim+1:3*h_dim,   h_dim+1:2*h_dim) =  F3
    Mmat( 2*h_dim+1:3*h_dim, 2*h_dim+1:3*h_dim) =  F0
    Mmat( 2*h_dim+1:3*h_dim, 3*h_dim+1:4*h_dim) = -F1
    Mmat( 3*h_dim+1:4*h_dim,           1:  h_dim) =  F3
    Mmat( 3*h_dim+1:4*h_dim,   h_dim+1:2*h_dim) = -F2
    Mmat( 3*h_dim+1:4*h_dim, 2*h_dim+1:3*h_dim) =  F1
    Mmat( 3*h_dim+1:4*h_dim, 3*h_dim+1:4*h_dim) =  F0

    ! -- 3. harmonic gradient at boundary nodes -----------------------
    deallocate(F, Fx, Fy, Fz)
    allocate(F(nbd, h_dim), Fx(nbd, h_dim), Fy(nbd, h_dim), Fz(nbd, h_dim))
    F = 0.0_r64; Fx = 0.0_r64; Fy = 0.0_r64; Fz = 0.0_r64
    call evaltensorproductharmonicgrad_r64(nbd, sxbd, order, Fx, Fy, Fz, F, ijidx)

    ! -- 4. q^{nm} kernel-coefficient builders (lptype 'd', 4 slots) --
    allocate(q_i(nbd, h_dim, 5), q_j(nbd, h_dim, 5), q_k(nbd, h_dim, 5))
    allocate(onm0(nbd, h_dim, 4), onm1(nbd, h_dim, 4), &
             onm2(nbd, h_dim, 4), onm3(nbd, h_dim, 4))

    allocate(dr(3, nbd))
    do k = 1, nbd
      dr(1, k) = stangbd(1, k) * swbd(k)
      dr(2, k) = stangbd(2, k) * swbd(k)
      dr(3, k) = stangbd(3, k) * swbd(k)
    end do

    call qak_qnm_i_r64(0_8, nbd, h_dim, sxbd, QAK_LPTYPE_D, &
                       F, Fx, Fy, Fz, dr, q_i, q_j, q_k)   ! gradxyz unused
    call qao_omeganm_i_r64(nbd, h_dim, 4_8, sxbd, dr,                    &
                           q_i(:,:,1:4), q_j(:,:,1:4), q_k(:,:,1:4), onm0)

    call qak_qnm_i_r64(1_8, nbd, h_dim, sxbd, QAK_LPTYPE_D, F, Fx, Fy, Fz, dr, q_i, q_j, q_k)
    call qao_omeganm_i_r64(nbd, h_dim, 4_8, sxbd, dr, q_i(:,:,1:4), q_j(:,:,1:4), q_k(:,:,1:4), onm1)

    call qak_qnm_i_r64(2_8, nbd, h_dim, sxbd, QAK_LPTYPE_D, F, Fx, Fy, Fz, dr, q_i, q_j, q_k)
    call qao_omeganm_i_r64(nbd, h_dim, 4_8, sxbd, dr, q_i(:,:,1:4), q_j(:,:,1:4), q_k(:,:,1:4), onm2)

    call qak_qnm_i_r64(3_8, nbd, h_dim, sxbd, QAK_LPTYPE_D, F, Fx, Fy, Fz, dr, q_i, q_j, q_k)
    call qao_omeganm_i_r64(nbd, h_dim, 4_8, sxbd, dr, q_i(:,:,1:4), q_j(:,:,1:4), q_k(:,:,1:4), onm3)

    ! -- 5. moments ----------------------------------------------------
    ! eval_moments_funvals_r64(..., order_arg, funvals_pre) writes
    !   funvals_pre(nbd, 2*(order_arg+1), m) = packed [N | M] moments.
    ! For Lap3dDLP_closepanel we need the M-block at moment count
    !   ncol = 2*(order+1).  Following lqs_eval_moments_funvals_mex:
    ! run the inner sub with moment_order = 2*order+1 (which gives
    ! 2*ncol total cols), then slice cols ncol+1 : 2*ncol (= M-block).
    allocate(M_all(nbd, ncol, m_tgt))
    M_all = 0.0_r64
    if (if_adapt) then
      block
        real(r64), allocatable :: funvals_full(:,:,:)
        integer(8) :: moment_order
        moment_order = 2_8*order + 1_8
        allocate(funvals_full(nbd, 2_8*ncol, m_tgt))
        call eval_moments_funvals_r64(m_tgt, t_x, nbd, sxbd, nquad, &
                                      moment_order, funvals_full)
        M_all = funvals_full(:, ncol+1_8:2_8*ncol, :)
        deallocate(funvals_full)
      end block
    else
      write(*,*) 'Lap3dDLP_closepanel_r64: if_adapt=.false. not yet implemented'
      write(*,*) '  (need a Fortran port of momentsallplain in moments.f).'
      error stop
    end if

    ! -- 6. assemble Omega_all -----------------------------------------
    allocate(Omega_all(m_tgt, 4*h_dim))
    call qao_omegaall_r64(m_tgt, nbd*ncol, nbd, h_dim, morder, t_x, &
                          reshape(M_all, (/nbd*ncol, m_tgt/)),       &
                          reshape(onm0, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm1, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm2, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm3, (/nbd*h_dim, 4_8/)),         &
                          ijidx, Omega_all)

    ! -- 7. solve Mmat^T * X^T = Omega_all^T, take first h_dim rows ----
    ! MATLAB does: A = (Mmat' \ Omega_all')' * [I; 0_{3h_dim}]
    !          = first h_dim columns of (Mmat' \ Omega_all')'
    !          = first h_dim rows transposed of (Mmat' \ Omega_all')
    ! r64 solve: build A_solve = Mmat^T, RHS = Omega_all^T, lu_solve...
    ! For r64 we *don't* use lu_solve (r64 in koorn_geom_mod is private);
    ! instead use LAPACK if available, or fall back to lu_solve_r128 path
    ! cast.  Since this module is library-only (no -framework Accelerate
    ! linkage assumed here), we ship a local copy of the partial-pivot
    ! LU below; same algorithm as lu_solve_r128 but at r64.
    allocate(rhs(4*h_dim, m_tgt))
    rhs = transpose(Omega_all)
    block
      real(r64) :: Asys(4*h_dim, 4*h_dim)
      Asys = transpose(Mmat)
      call lu_solve_local_r64(4_8*h_dim, Asys, m_tgt, rhs)
    end block
    ! Ac(m_tgt, npat) = rhs(1:h_dim, :)^T   (only the first h_dim block)
    do k = 1, m_tgt
      do kk = 1, npat
        Ac(k, kk) = rhs(kk, k)
      end do
    end do

  end subroutine Lap3dDLP_closepanel_r64

  ! ================================================================
  ! Lap3dDLP_closepanel_r128
  !
  ! r128 sibling.  Open r128 dependencies (NOT YET PORTED):
  !   - evaltensorproductharmonicgrad_r128  (not in QA-legacy)
  !   - eval_moments_funvals_r128            (not in LineQuaaadrature-legacy)
  ! Both currently error_stop.  Everything else along the chain is
  ! already available at r128 in the legacy packages and is wired here.
  ! ================================================================
  subroutine Lap3dDLP_closepanel_r128(m_tgt, t_x, npat, s_x, order, ref, &
                                      if_adapt, Ac)
    integer(8), intent(in)  :: m_tgt, npat, order, ref
    real(r128), intent(in)  :: t_x(3, m_tgt)
    real(r128), intent(in)  :: s_x(3, npat)
    logical,    intent(in)  :: if_adapt
    real(r128), intent(out) :: Ac(m_tgt, npat)

    ! Silence unused-arg warnings while the r128 chain has gaps.
    Ac = 0.0_r128
    if (.false.) then
      write(*,*) m_tgt, npat, order, ref, if_adapt, t_x(1,1), s_x(1,1)
    end if

    write(*,*) 'Lap3dDLP_closepanel_r128: not yet runnable.'
    write(*,*) '  Missing r128 ports:'
    write(*,*) '    1) evaltensorproductharmonicgrad_r128 (QA-legacy harmonic_mod)'
    write(*,*) '    2) eval_moments_funvals_r128         (LineQuaaadrature-legacy solidangle_mod)'
    write(*,*) '  Once those land, mirror Lap3dDLP_closepanel_r64 with kind(r128).'
    error stop
  end subroutine Lap3dDLP_closepanel_r128

  ! ================================================================
  ! lu_solve_local_r64  (private)
  ! In-place LU with partial pivoting.  Same algorithm as
  ! koorn_geom_mod::lu_solve, but local to avoid making the QA-legacy
  ! r64 lu_solve public.  Sized for the small block matrices (~64x64)
  ! used here.
  ! ================================================================
  subroutine lu_solve_local_r64(n, A, k, B)
    integer(8), intent(in)    :: n, k
    real(r64),  intent(inout) :: A(n,n), B(n,k)
    integer(8) :: i, j, p, col, pivot_row
    real(r64)  :: pivot, tmp, factor, rowA(n), rowB(k)

    do col = 1, n
      pivot_row = col
      pivot     = abs(A(col,col))
      do i = col+1, n
        if (abs(A(i,col)) > pivot) then
          pivot     = abs(A(i,col))
          pivot_row = i
        end if
      end do
      if (pivot_row /= col) then
        rowA = A(col,:);  A(col,:) = A(pivot_row,:);  A(pivot_row,:) = rowA
        rowB = B(col,:);  B(col,:) = B(pivot_row,:);  B(pivot_row,:) = rowB
      end if
      do i = col+1, n
        factor   = A(i,col) / A(col,col)
        A(i,col) = factor
        do j = col+1, n
          A(i,j) = A(i,j) - factor*A(col,j)
        end do
        do p = 1, k
          B(i,p) = B(i,p) - factor*B(col,p)
        end do
      end do
    end do
    do col = n, 1, -1
      do p = 1, k
        tmp = B(col,p)
        do j = col+1, n
          tmp = tmp - A(col,j)*B(j,p)
        end do
        B(col,p) = tmp / A(col,col)
      end do
    end do
  end subroutine lu_solve_local_r64

end module lap3d_close_mod
