#ifdef BIESOLVER_WASM_SCALAR_ONLY
#define cpu_time biesolver_scalar_noop_cpu_time
#endif
module lap3d_close_mod

  use quatapproximation_mod, only: r64, r128, gauss_r64
#ifdef BIESOLVER_STELLARATOR_BUILD
  use iso_c_binding, only: c_char
#endif
#ifndef BIESOLVER_R64_ONLY
  use quatapproximation_mod, only: gauss_r128
#endif
  use harmonic_mod, only: evaltensorproductharmonicgrad_r64
  use tensor_geom_mod, only: line3quadr_3dline_T
#ifndef BIESOLVER_R64_ONLY
  use tensor_geom_mod, only: line3quadr_3dline_T_r128
#endif
  use qkernel_mod, only: qak_qnm_i_r64, QAK_LPTYPE_D
#ifndef BIESOLVER_R64_ONLY
  use qkernel_mod, only: qak_qnm_i_r128
#endif
  use omega_mod, only: qao_omeganm_i_r64, qao_omegaall_r64
#ifndef BIESOLVER_R64_ONLY
  use omega_mod, only: qao_omeganm_i_r128, qao_omegaall_r128
#endif
#ifndef BIESOLVER_R64_ONLY
  use koorn_geom_mod, only: lu_solve_r128
#endif
  use solidangle_mod, only: eval_moments_funvals_r64

  implicit none
  private
  public :: Lap3dDLP_closepanel_r64
#ifndef BIESOLVER_R64_ONLY
  public :: Lap3dDLP_closepanel_r128
#endif
  public :: build_tc_chialpha_r64
  public :: build_closepanel_precomp_r64
  public :: rrq_r64

  integer(8), parameter :: SBDNP = 8_8

#ifndef BIESOLVER_WASM_SCALAR_ONLY
  integer(8), save :: CREF_nquad = -1_8, CREF_korder = -1_8
  integer(8), save :: CREF_kpols = -1_8
  logical, save :: CREF_set = .false.
  real(r64), allocatable, save :: CREF_tgl(:), CREF_wgl(:), CREF_Dgl(:,:)
  real(r64), allocatable, save :: CREF_w_bclag(:), CREF_Legmat(:,:)
  real(r64), allocatable, save :: CREF_umatr(:,:), CREF_vmatr(:,:)
#endif

  interface
#ifdef BIESOLVER_STELLARATOR_BUILD
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc) bind(C, name="biesolver_dgemm")
      import r64, c_char
      character(kind=c_char), value :: transa, transb
      integer(8), intent(in) :: m, n, k, lda, ldb, ldc
      real(r64), intent(in) :: alpha, beta
      real(r64), intent(in) :: a(lda,*), b(ldb,*)
      real(r64), intent(inout) :: c(ldc,*)
    end subroutine dgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc) bind(C, name="biesolver_zgemm")
      import r64, c_char
      character(kind=c_char), value :: transa, transb
      integer(8), intent(in) :: m, n, k, lda, ldb, ldc
      complex(r64), intent(in) :: alpha, beta
      complex(r64), intent(in) :: a(lda,*), b(ldb,*)
      complex(r64), intent(inout) :: c(ldc,*)
    end subroutine zgemm

    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info) &
        bind(C, name="biesolver_dgesv")
      import r64
      integer(8), intent(in) :: n, nrhs, lda, ldb
      real(r64), intent(inout) :: a(lda,*), b(ldb,*)
      integer(8), intent(out) :: ipiv(*)
      integer(8), intent(out) :: info
    end subroutine dgesv
#else
    subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc)
      import r64
      character(len=1), intent(in) :: transa, transb
      integer(8), intent(in) :: m, n, k, lda, ldb, ldc
      real(r64), intent(in) :: alpha, beta
      real(r64), intent(in) :: a(lda,*), b(ldb,*)
      real(r64), intent(inout) :: c(ldc,*)
    end subroutine dgemm

    subroutine zgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
                     beta, c, ldc)
      import r64
      character(len=1), intent(in) :: transa, transb
      integer(8), intent(in) :: m, n, k, lda, ldb, ldc
      complex(r64), intent(in) :: alpha, beta
      complex(r64), intent(in) :: a(lda,*), b(ldb,*)
      complex(r64), intent(inout) :: c(ldc,*)
    end subroutine zgemm

    subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
      import r64
      integer(8), intent(in) :: n, nrhs, lda, ldb
      real(r64), intent(inout) :: a(lda,*), b(ldb,*)
      integer(8), intent(out) :: ipiv(*)
      integer(8), intent(out) :: info
    end subroutine dgesv
#endif
  end interface

contains

#ifdef BIESOLVER_WASM_SCALAR_ONLY
  subroutine biesolver_scalar_noop_cpu_time(value)
    real(r64), intent(out) :: value
    value = 0.0_r64
  end subroutine biesolver_scalar_noop_cpu_time
#endif

  subroutine close_ref_build_r64(nquad, korder, kpols, tgl, wgl, Dgl, &
                                 w_bclag, Legmat, umatr, vmatr)
    use koorn_geom_mod,        only: koorn_vals2coefs_coefs2vals
    use quatapproximation_mod, only: gauss_r64, bclaginterpweights_r64
    use linequaaadrature_mod,  only: legeexps_r64
    integer(8), intent(in) :: nquad, korder, kpols
    real(r64), intent(out) :: tgl(nquad), wgl(nquad)
    real(r64), intent(out) :: Dgl(nquad,nquad), w_bclag(nquad)
    real(r64), intent(out) :: Legmat(nquad,nquad)
    real(r64), intent(out) :: umatr(kpols,kpols), vmatr(kpols,kpols)
    real(r64) :: vtmp(nquad,nquad)

    call gauss_r64(nquad, tgl, wgl, Dgl)
    call bclaginterpweights_r64(nquad, tgl, w_bclag)
    call legeexps_r64(2_8, nquad, tgl, Legmat, vtmp, wgl)

    umatr = 0.0_r64;  vmatr = 0.0_r64
    call koorn_vals2coefs_coefs2vals(korder, kpols, umatr, vmatr)
  end subroutine close_ref_build_r64

#ifndef BIESOLVER_WASM_SCALAR_ONLY
  subroutine close_ref_ensure_r64(nquad, korder, kpols)
    integer(8), intent(in) :: nquad, korder, kpols

    if (CREF_set .and. CREF_nquad == nquad .and. &
        CREF_korder == korder .and. CREF_kpols == kpols) return

    !$omp critical (close_ref_build)
    if (.not. (CREF_set .and. CREF_nquad == nquad .and. &
               CREF_korder == korder .and. CREF_kpols == kpols)) then

    if (CREF_set) deallocate(CREF_tgl, CREF_wgl, CREF_Dgl, &
                             CREF_w_bclag, CREF_Legmat, &
                             CREF_umatr, CREF_vmatr)
    allocate(CREF_tgl(nquad), CREF_wgl(nquad), CREF_Dgl(nquad,nquad))
    allocate(CREF_w_bclag(nquad), CREF_Legmat(nquad,nquad))
    allocate(CREF_umatr(kpols,kpols), CREF_vmatr(kpols,kpols))
    call close_ref_build_r64(nquad, korder, kpols, CREF_tgl, CREF_wgl, &
                             CREF_Dgl, CREF_w_bclag, CREF_Legmat, &
                             CREF_umatr, CREF_vmatr)
    CREF_nquad = nquad;  CREF_korder = korder;  CREF_kpols = kpols
    CREF_set = .true.

    end if
    !$omp end critical (close_ref_build)

  end subroutine close_ref_ensure_r64
#endif

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

    allocate(F(npat, h_dim), Fx(npat, h_dim), Fy(npat, h_dim), Fz(npat, h_dim))
    allocate(ijidx(2, h_dim))
    F = 0.0_r64; Fx = 0.0_r64; Fy = 0.0_r64; Fz = 0.0_r64; ijidx = 0_8
    call evaltensorproductharmonicgrad_r64(npat, s_x, order, Fx, Fy, Fz, F, ijidx)

    allocate(F0(h_dim, h_dim), F1(h_dim, h_dim), F2(h_dim, h_dim), F3(h_dim, h_dim))
    F0 = 0.0_r64
    F1 = Fx ; F2 = Fy ; F3 = Fz   ! 1st h_dim source-nodes (assumes npat=h_dim)

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

    deallocate(F, Fx, Fy, Fz)
    allocate(F(nbd, h_dim), Fx(nbd, h_dim), Fy(nbd, h_dim), Fz(nbd, h_dim))
    F = 0.0_r64; Fx = 0.0_r64; Fy = 0.0_r64; Fz = 0.0_r64
    call evaltensorproductharmonicgrad_r64(nbd, sxbd, order, Fx, Fy, Fz, F, ijidx)

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

    allocate(Omega_all(m_tgt, 4*h_dim))
    call qao_omegaall_r64(m_tgt, nbd*ncol, nbd, h_dim, morder, t_x, &
                          reshape(M_all, (/nbd*ncol, m_tgt/)),       &
                          reshape(onm0, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm1, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm2, (/nbd*h_dim, 4_8/)),         &
                          reshape(onm3, (/nbd*h_dim, 4_8/)),         &
                          ijidx, Omega_all)

    allocate(rhs(4*h_dim, m_tgt))
    rhs = transpose(Omega_all)
    block
      real(r64) :: Asys(4*h_dim, 4*h_dim)
      Asys = transpose(Mmat)
      call lu_solve_local_r64(4_8*h_dim, Asys, m_tgt, rhs)
    end block
    do k = 1, m_tgt
      do kk = 1, npat
        Ac(k, kk) = rhs(kk, k)
      end do
    end do

  end subroutine Lap3dDLP_closepanel_r64

#ifndef BIESOLVER_R64_ONLY
  subroutine Lap3dDLP_closepanel_r128(m_tgt, t_x, npat, s_x, order, ref, &
                                      if_adapt, Ac)
    integer(8), intent(in)  :: m_tgt, npat, order, ref
    real(r128), intent(in)  :: t_x(3, m_tgt)
    real(r128), intent(in)  :: s_x(3, npat)
    logical,    intent(in)  :: if_adapt
    real(r128), intent(out) :: Ac(m_tgt, npat)

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
#endif

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

  subroutine build_tc_chialpha_r64(m, r0, nbd, sbdnp, nquad, ncoeff, &
                                   sxbd, stangbd, sspbd, w1, w3, &
                                   Fx, Fy, Fz, IalphaAsvestas, Ichi, Ialpha)
    integer(8),   intent(in)  :: m, nbd, sbdnp, nquad, ncoeff
    real(r64),    intent(in)  :: r0(3,m)
    real(r64),    intent(in)  :: sxbd(3,nbd), stangbd(3,nbd), sspbd(nbd)
    real(r64),    intent(in)  :: w1(nquad,sbdnp,m), w3(nquad,sbdnp,m)
    complex(r64), intent(in)  :: Fx(nbd,ncoeff), Fy(nbd,ncoeff), Fz(nbd,ncoeff)
    real(r64),    intent(in)  :: IalphaAsvestas(m)
    complex(r64), intent(out) :: Ichi(m,ncoeff,4), Ialpha(m,ncoeff,4)

    complex(r64), parameter :: ONEC = ( 1.0_r64, 0.0_r64)
    complex(r64), parameter :: IMA  = ( 0.0_r64, 1.0_r64)
    real(r64),    parameter :: ONE  =  1.0_r64
    real(r64),    parameter :: MONE = -1.0_r64

    real(r64) :: dxdt(nbd), dydt(nbd), dzdt(nbd)
    real(r64) :: rx(nbd,m), ry(nbd,m), rz(nbd,m)
    real(r64) :: xx(nbd,m), yy(nbd,m), zz(nbd,m)
    real(r64) :: xy(nbd,m), xz(nbd,m), yz(nbd,m)
    real(r64) :: rinvw1(nbd,m), r3invw3(nbd,m), rr(nbd,m)
    real(r64) :: ax(nbd,m), ay(nbd,m), az(nbd,m)
    real(r64) :: bx(nbd,m), by(nbd,m), bz3(nbd,m)
    real(r64) :: g1(nbd,m), g2(nbd,m), g3(nbd,m)
    real(r64) :: Fxri(nbd,2*ncoeff), Fyri(nbd,2*ncoeff)
    real(r64) :: Fzri(nbd,2*ncoeff)
    real(r64) :: Ichi_ri(m,2*ncoeff), Ialpha_ri(m,2*ncoeff,4)
    real(r64) :: c0(m), c1(m), c2(m), c3(m)
    complex(r64) :: Fx1(2), Fy1(2), Fz1(2)
    integer(8) :: i, j, k, ell, iq

    dxdt = stangbd(1,:)*sspbd
    dydt = stangbd(2,:)*sspbd
    dzdt = stangbd(3,:)*sspbd

    do j = 1, m
      rx(:,j) = sxbd(1,:) - r0(1,j)
      ry(:,j) = sxbd(2,:) - r0(2,j)
      rz(:,j) = sxbd(3,:) - r0(3,j)
    end do
    xx = rx*rx;  yy = ry*ry;  zz = rz*rz
    xy = rx*ry;  xz = rx*rz;  yz = ry*rz

    rr = 1.0_r64/sqrt(xx + yy + zz)
    do j = 1, m
      do ell = 1, sbdnp
        do iq = 1, nquad
          i = (ell-1_8)*nquad + iq
          rinvw1(i,j)  = rr(i,j)*w1(iq,ell,j)
          r3invw3(i,j) = rr(i,j)*rr(i,j)*rr(i,j)*w3(iq,ell,j)
        end do
      end do
    end do

    do j = 1, m
      ax(:,j)  = dxdt*rinvw1(:,j)
      ay(:,j)  = dydt*rinvw1(:,j)
      az(:,j)  = dzdt*rinvw1(:,j)
      bx(:,j)  = dxdt*r3invw3(:,j)
      by(:,j)  = dydt*r3invw3(:,j)
      bz3(:,j) = dzdt*r3invw3(:,j)
    end do

    Fxri(:,1:ncoeff) = real(Fx, r64)
    Fxri(:,ncoeff+1:2*ncoeff) = aimag(Fx)
    Fyri(:,1:ncoeff) = real(Fy, r64)
    Fyri(:,ncoeff+1:2*ncoeff) = aimag(Fy)
    Fzri(:,1:ncoeff) = real(Fz, r64)
    Fzri(:,ncoeff+1:2*ncoeff) = aimag(Fz)

    Ichi = (0.0_r64, 0.0_r64)
    Ichi_ri = 0.0_r64
    g1 = -rz*ay + ry*az
    g2 =  rz*ax - rx*az
    g3 = -ry*ax + rx*ay

    call dgemm('T','N', m, 2*ncoeff, nbd, ONE, g1, nbd, Fxri, nbd, &
               ONE, Ichi_ri, m)
    call dgemm('T','N', m, 2*ncoeff, nbd, ONE, g2, nbd, Fyri, nbd, &
               ONE, Ichi_ri, m)
    call dgemm('T','N', m, 2*ncoeff, nbd, ONE, g3, nbd, Fzri, nbd, &
               ONE, Ichi_ri, m)
    Ichi(:,:,1) = cmplx(Ichi_ri(:,1:ncoeff), &
                        Ichi_ri(:,ncoeff+1:2*ncoeff), r64)
    Ichi(:,1,1) = cmplx(IalphaAsvestas, 0.0_r64, r64)

    Ialpha = (0.0_r64, 0.0_r64)
    Ialpha_ri = 0.0_r64
    g1 = -(yy + zz)*bx + xy*by + xz*bz3
    g2 =  xy*bx - (xx + zz)*by + yz*bz3
    g3 =  xz*bx + yz*by - (xx + yy)*bz3

    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g1, nbd, Fxri, nbd, &
               ONE, Ialpha_ri(:,:,1), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g1, nbd, Fzri, nbd, &
               ONE, Ialpha_ri(:,:,3), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, MONE, g1, nbd, Fyri, nbd, &
               ONE, Ialpha_ri(:,:,4), m)

    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g2, nbd, Fyri, nbd, &
               ONE, Ialpha_ri(:,:,1), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, MONE, g2, nbd, Fzri, nbd, &
               ONE, Ialpha_ri(:,:,2), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g2, nbd, Fxri, nbd, &
               ONE, Ialpha_ri(:,:,4), m)

    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g3, nbd, Fzri, nbd, &
               ONE, Ialpha_ri(:,:,1), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, ONE,  g3, nbd, Fyri, nbd, &
               ONE, Ialpha_ri(:,:,2), m)
    call dgemm('T','N', m, 2*ncoeff, nbd, MONE, g3, nbd, Fxri, nbd, &
               ONE, Ialpha_ri(:,:,3), m)

    do k = 1, 4
      Ialpha(:,:,k) = cmplx(Ialpha_ri(:,1:ncoeff,k), &
                            Ialpha_ri(:,ncoeff+1:2*ncoeff,k), r64)
    end do

    Fx1 = [ (0.0_r64,0.0_r64), ONEC ]
    Fy1 = [ cmplx(-sqrt(2.0_r64)/2.0_r64, 0.0_r64, r64), (0.0_r64,0.0_r64) ]
    Fz1 = [ sqrt(2.0_r64)/2.0_r64*IMA, (0.0_r64,0.0_r64) ]

    c0 = IalphaAsvestas
    c1 = matmul(transpose(rinvw1), dxdt)
    c2 = matmul(transpose(rinvw1), dydt)
    c3 = matmul(transpose(rinvw1), dzdt)

    do j = 1, 2
      Ialpha(:,j+1,1) = c1*Fx1(j) + c2*Fy1(j) + c3*Fz1(j)
      Ialpha(:,j+1,2) = c0*Fx1(j) - c2*Fz1(j) + c3*Fy1(j)
      Ialpha(:,j+1,3) = c0*Fy1(j) + c1*Fz1(j) - c3*Fx1(j)
      Ialpha(:,j+1,4) = c0*Fz1(j) - c1*Fy1(j) + c2*Fx1(j)
    end do

  end subroutine build_tc_chialpha_r64

  subroutine build_closepanel_precomp_r64(n, sx, snx, sw, r_vert, &
                                          nterms, h_dim, nbd, sbdnp, nquad, &
                                          alpha, exterior, ichart, sxbd_chart, &
                                          tgl, wgl, Dgl, w_bclag, Legmat, &
                                          R, c, sxbd, swbd, stangbd, sspbd, &
                                          qhat, Fbd, Fxbd, Fybd, Fzbd, Mmatrix)
    use koorn_geom_mod, only: circumcircle_transform_3d, &
                              koorn_vals2coefs_coefs2vals, line3quadr_3dline
    use harmonic_mod,   only: l3dtavecevalmat_r64
    use quatapproximation_mod, only: gauss_r64, bclaginterpweights_r64
    use linequaaadrature_mod,  only: legeexps_r64
    integer(8), intent(in)  :: n, nterms, h_dim, nbd, sbdnp, nquad
    real(r64),  intent(in)  :: sx(3,n), snx(3,n), sw(n), r_vert(3,3), alpha
    integer(8), intent(in)  :: ichart
    real(r64),  intent(in)  :: sxbd_chart(3,nbd)
    logical,    intent(in)  :: exterior
    real(r64),  intent(out) :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad)
    real(r64),  intent(out) :: w_bclag(nquad), Legmat(nquad,nquad)
    real(r64),  intent(out) :: R(3,3), c(3)
    real(r64),  intent(out) :: sxbd(3,nbd), swbd(nbd)
    real(r64),  intent(out) :: stangbd(3,nbd), sspbd(nbd)
    real(r64),  intent(out) :: qhat(3)
    complex(r64), intent(out) :: Fbd(nbd,(nterms+1)**2), Fxbd(nbd,(nterms+1)**2)
    complex(r64), intent(out) :: Fybd(nbd,(nterms+1)**2), Fzbd(nbd,(nterms+1)**2)
    real(r64),  intent(out) :: Mmatrix(4*h_dim,4*h_dim)

    real(r64), parameter :: OFFSET_FACTOR = 1.25_r64   ! rrq linequadv2.f:8140

    real(r64)  :: vtmp(nquad,nquad), umatr(n,n), vmatr(n,n)
    real(r64)  :: c0(3), alpha_circ, sxt(3,n), snxt(3,n)
    real(r64)  :: tpan(sbdnp+1), r_vert_local(3,3), sx3min, pi
    real(r64)  :: F0(n,h_dim), F1(n,h_dim), F2(n,h_dim), F3(n,h_dim)
    real(r64)  :: MmatrixS(h_dim,h_dim)
    integer(8) :: korder, kpols, idxvec(h_dim), i, k, ij, t1, t2, ier

    pi = 4.0_r64*atan(1.0_r64)

    korder = nterms - 1_8
    kpols  = nterms*(nterms+1_8)/2_8
#ifdef BIESOLVER_WASM_SCALAR_ONLY
    call close_ref_build_r64(nquad, korder, kpols, tgl, wgl, Dgl, &
                             w_bclag, Legmat, umatr, vmatr)
#else
    call close_ref_ensure_r64(nquad, korder, kpols)
    tgl = CREF_tgl;  wgl = CREF_wgl;  Dgl = CREF_Dgl
    w_bclag = CREF_w_bclag;  Legmat = CREF_Legmat
    umatr = CREF_umatr;  vmatr = CREF_vmatr
#endif

    call circumcircle_transform_3d(r_vert, R, c0, alpha_circ)
    do i = 1, n
      sxt(:,i) = alpha*matmul(R, sx(:,i) - c0)
    end do
    sx3min = max(0.0_r64, -minval(sxt(3,:)))
    c      = c0 - (OFFSET_FACTOR*sx3min/alpha)*R(3,:)
    do i = 1, n
      sxt(:,i)  = alpha*matmul(R, sx(:,i) - c)
      snxt(:,i) = matmul(R, snx(:,i))
    end do

    do k = 1, sbdnp+1
      tpan(k) = real(k-1, r64)*2.0_r64*pi/real(sbdnp, r64)
    end do
    sxbd = 0.0_r64;  swbd = 0.0_r64;  stangbd = 0.0_r64;  sspbd = 0.0_r64
    r_vert_local = 0.0_r64
    call line3quadr_3dline(sxt, korder, kpols, umatr, nquad, tgl, wgl, Dgl, &
                           sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, &
                           r_vert_local)

    if (ichart == 1_8) then
      block
        real(r64) :: sxp(3,nbd), pt, sw0(nbd)
        integer(8) :: ell, i0, i1, iq
        do i = 1, nbd
          sxbd(:,i) = alpha*matmul(R, sxbd_chart(:,i) - c)
        end do
        pt = 2.0_r64*pi/real(sbdnp, r64)
        do ell = 1, sbdnp
          i0 = (ell-1_8)*nquad + 1_8
          i1 = ell*nquad
          sxp(1,i0:i1) = (2.0_r64/pt)*matmul(Dgl, sxbd(1,i0:i1))
          sxp(2,i0:i1) = (2.0_r64/pt)*matmul(Dgl, sxbd(2,i0:i1))
          sxp(3,i0:i1) = (2.0_r64/pt)*matmul(Dgl, sxbd(3,i0:i1))
          do iq = 1, nquad
            sw0(i0+iq-1_8) = 0.5_r64*wgl(iq)*pt
          end do
        end do
        sspbd = sqrt(sxp(1,:)**2 + sxp(2,:)**2 + sxp(3,:)**2)
        stangbd(1,:) = sxp(1,:)/sspbd
        stangbd(2,:) = sxp(2,:)/sspbd
        stangbd(3,:) = sxp(3,:)/sspbd
        swbd = sw0*sspbd
      end block
    end if

    sspbd = (pi/real(sbdnp, r64))*sspbd

    qhat(1) = sum(snxt(1,:))/real(n, r64)
    qhat(2) = sum(snxt(2,:))/real(n, r64)
    qhat(3) = sum(snxt(3,:))/real(n, r64)
    if (.not. exterior) qhat = -qhat
    qhat = qhat/sqrt(dot_product(qhat, qhat))

    block
      real(r64)    :: p(3,nbd)
      complex(r64) :: Gx(nbd,(nterms+1)**2), Gy(nbd,(nterms+1)**2)
      complex(r64) :: Gz(nbd,(nterms+1)**2)
      p(1,:) = sxbd(2,:);  p(2,:) = sxbd(3,:);  p(3,:) = sxbd(1,:)
      ier = 0
      call l3dtavecevalmat_r64(p, nbd, nterms, Fbd, Gx, Gy, Gz, ier)
      Fxbd = Gz;  Fybd = Gx;  Fzbd = Gy
    end block

    t1 = 0_8;  t2 = 0_8
    do ij = 0, nterms
      do k = -ij, ij
        t2 = t2 + 1_8
        if (ij > 0_8 .and. k > 0_8) then
          t1 = t1 + 1_8;  idxvec(t1) = t2
        end if
      end do
    end do

    block
      real(r64)    :: p(3,n)
      complex(r64) :: Fc(n,(nterms+1)**2), Gx(n,(nterms+1)**2)
      complex(r64) :: Gy(n,(nterms+1)**2), Gz(n,(nterms+1)**2)
      p(1,:) = sxt(2,:);  p(2,:) = sxt(3,:);  p(3,:) = sxt(1,:)
      ier = 0
      call l3dtavecevalmat_r64(p, n, nterms, Fc, Gx, Gy, Gz, ier)
      do k = 1, h_dim
        F0(:,k) = aimag(Fc(:,idxvec(k)))
        F1(:,k) = aimag(Gz(:,idxvec(k)))
        F2(:,k) = aimag(Gx(:,idxvec(k)))
        F3(:,k) = aimag(Gy(:,idxvec(k)))
      end do
    end block

    Mmatrix = 0.0_r64
    Mmatrix(     1:  n, (  h_dim+1):(2*h_dim)) = -F1
    Mmatrix(     1:  n, (2*h_dim+1):(3*h_dim)) = -F2
    Mmatrix(     1:  n, (3*h_dim+1):(4*h_dim)) = -F3
    Mmatrix(  n+1:2*n,         1   :   h_dim ) =  F1
    Mmatrix(  n+1:2*n, (2*h_dim+1):(3*h_dim)) = -F3
    Mmatrix(  n+1:2*n, (3*h_dim+1):(4*h_dim)) =  F2
    Mmatrix(2*n+1:3*n,         1   :   h_dim ) =  F2
    Mmatrix(2*n+1:3*n, (  h_dim+1):(2*h_dim)) =  F3
    Mmatrix(2*n+1:3*n, (3*h_dim+1):(4*h_dim)) = -F1
    Mmatrix(3*n+1:4*n,         1   :   h_dim ) =  F3
    Mmatrix(3*n+1:4*n, (  h_dim+1):(2*h_dim)) = -F2
    Mmatrix(3*n+1:4*n, (2*h_dim+1):(3*h_dim)) =  F1

    do k = 1, h_dim
      MmatrixS(:,k) = snxt(1,:)*F1(:,k) + snxt(2,:)*F2(:,k) + snxt(3,:)*F3(:,k)
    end do
    Mmatrix(      1:  h_dim,       1:  h_dim) = MmatrixS
    Mmatrix(h_dim+1:2*h_dim, h_dim+1:2*h_dim) = F0

  end subroutine build_closepanel_precomp_r64

  subroutine rrq_r64(m, tx, n, sx, snx, sw, rts, rps, &
                     order, nquad, orderff, distff, exterior, isimd, &
                     ichart, sxbd_chart, rv_chart, &
                     As, Ad, Omega, IalphaAsvestas, timeinfo)
    use patch_refine_mod, only: patch_levels_t, build_patch_levels_r64, &
                                clear_patch_levels_r64, &
                                lap3dsdlpmat_levels_r64
    use koorn_geom_mod,   only: koorn_vals2coefs_coefs2vals, line3quadr_3dline
    use lq_kernel_mod,    only: build_ssq_weights_r64
    use solidangle_mod,   only: evaluate_solid_angle_integral_fast_r64
    use omega_mod,        only: qao_omegasdlp_r64
    integer(8), intent(in)    :: m, n, order, nquad, orderff, isimd
    real(r64),  intent(in)    :: tx(3,m), sx(3,n), snx(3,n), sw(n)
    real(r64),  intent(in)    :: rts(3,n), rps(3,n), distff
    logical,    intent(in)    :: exterior
    real(r64),  intent(out)   :: As(m,n), Ad(m,n)
    real(r64),  intent(out)   :: Omega(4*(order*(order+1)/2),m)
    real(r64),  intent(out)   :: IalphaAsvestas(m)
    real(r64),  intent(inout) :: timeinfo(20)
    integer(8), intent(in)    :: ichart
    real(r64),  intent(in)    :: sxbd_chart(3,3*nquad), rv_chart(3,3)

    integer(8), parameter :: LEN1 = 4_8, LEN2 = 8_8, LEN3 = 16_8
    integer(8), parameter :: NLEVEL = 4_8
    real(r64),  parameter :: RHO_SSQ = 100.0_r64                 ! rrq :609

    integer(8) :: nterms, h_dim, ncoeff, sbdnp, nbd, korder, kpols
    integer(8) :: nb1, nb2, nb3, i, j, k, ij, t1, t2, ms
    real(r64)  :: alpha, pi, fac, t0, t1c
    type(patch_levels_t) :: lv

    real(r64) :: tgl(nquad), wgl(nquad), Dgl(nquad,nquad), w_bclag(nquad)
    real(r64) :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
    real(r64) :: R(3,3), c(3), qhat(3), r_vert(3,3)
    real(r64) :: sxbd(3,3*nquad), swbd(3*nquad)
    real(r64) :: sxbd_raw(3,3*nquad)
    real(r64) :: stangbd(3,3*nquad), sspbd(3*nquad)
    real(r64) :: sxbd1(3,3*LEN1*nquad), stangbd1(3,3*LEN1*nquad)
    real(r64) :: swbd1(3*LEN1*nquad)
    real(r64) :: sxbd2(3,3*LEN2*nquad), stangbd2(3,3*LEN2*nquad)
    real(r64) :: swbd2(3*LEN2*nquad)
    real(r64) :: sxbd3(3,3*LEN3*nquad), stangbd3(3,3*LEN3*nquad)
    real(r64) :: swbd3(3*LEN3*nquad)
    complex(r64) :: Fbd(3*nquad,(order+1)**2), Fxbd(3*nquad,(order+1)**2)
    complex(r64) :: Fybd(3*nquad,(order+1)**2), Fzbd(3*nquad,(order+1)**2)
    real(r64) :: Mmatrix(4*(order*(order+1)/2),4*(order*(order+1)/2))
    real(r64) :: sxt(3,n), dl(nquad), dr(nquad)
    real(r64) :: umatr(order*(order+1)/2,order*(order+1)/2)
    real(r64) :: vmatr(order*(order+1)/2,order*(order+1)/2)
    integer(8) :: idxs(m)

    pi     = 4.0_r64*atan(1.0_r64)
    nterms = order
    h_dim  = nterms*(nterms+1_8)/2_8
    ncoeff = (nterms+1_8)*(nterms+2_8)/2_8
    sbdnp  = 3_8                       ! rrq: sbdnp = 3*len0, len0 = 1
    nbd    = sbdnp*nquad
    korder = nterms - 1_8
    kpols  = nterms*(nterms+1_8)/2_8
    nb1    = 3_8*LEN1*nquad
    nb2    = 3_8*LEN2*nquad
    nb3    = 3_8*LEN3*nquad

    As = 0.0_r64;  Ad = 0.0_r64
    Omega = 0.0_r64;  IalphaAsvestas = 0.0_r64

    alpha = 1.25_r64/sqrt(2.0_r64*sum(sw))

    call cpu_time(t0)
    call build_patch_levels_r64(korder, n, sx, rts, rps, orderff, distff, &
                                NLEVEL, lv)
    call cpu_time(t1c);  timeinfo(1) = timeinfo(1) + (t1c - t0)

    call cpu_time(t0)

#ifdef BIESOLVER_WASM_SCALAR_ONLY
    call close_ref_build_r64(nquad, korder, kpols, tgl, wgl, Dgl, &
                             w_bclag, Legmat, umatr, vmatr)
#else
    call close_ref_ensure_r64(nquad, korder, kpols)
    tgl = CREF_tgl;  wgl = CREF_wgl;  Dgl = CREF_Dgl
    w_bclag = CREF_w_bclag;  Legmat = CREF_Legmat
    umatr = CREF_umatr;  vmatr = CREF_vmatr
#endif
    block
      real(r64) :: tp(sbdnp+1), bw(3*nquad)
      real(r64) :: bt(3,3*nquad), bs(3*nquad)
      do k = 1, sbdnp+1_8
        tp(k) = real(k-1, r64)*2.0_r64*pi/real(sbdnp, r64)
      end do
      sxbd_raw = 0.0_r64; bw = 0.0_r64; bt = 0.0_r64; bs = 0.0_r64
      r_vert = 0.0_r64
      call line3quadr_3dline(sx, korder, kpols, umatr, nquad, tgl, wgl, Dgl, &
                             sbdnp, tp, nbd, sxbd_raw, bw, bt, bs, r_vert)
    end block
    if (ichart == 1_8) sxbd_raw = sxbd_chart

    if (ichart == 1_8) r_vert = rv_chart
    call build_closepanel_precomp_r64(n, sx, snx, sw, r_vert, &
                                      nterms, h_dim, nbd, sbdnp, nquad, &
                                      alpha, exterior, ichart, sxbd_chart, &
                                      tgl, wgl, Dgl, w_bclag, Legmat, &
                                      R, c, sxbd, swbd, stangbd, sspbd, &
                                      qhat, Fbd, Fxbd, Fybd, Fzbd, Mmatrix)

    do i = 1, n
      sxt(:,i) = alpha*matmul(R, sx(:,i) - c)
    end do

    block
      real(r64) :: tp(3*LEN1+1), ss(3*LEN1*nquad), rv(3,3)
      do k = 1, 3_8*LEN1+1_8
        tp(k) = real(k-1, r64)*2.0_r64*pi/real(3_8*LEN1, r64)
      end do
      sxbd1 = 0.0_r64; swbd1 = 0.0_r64; stangbd1 = 0.0_r64; ss = 0.0_r64
      rv = 0.0_r64
      call line3quadr_3dline(sxt, korder, kpols, umatr, nquad, tgl, wgl, Dgl, &
                             3_8*LEN1, tp, nb1, sxbd1, swbd1, stangbd1, ss, rv)
    end block
    block
      real(r64) :: tp(3*LEN2+1), ss(3*LEN2*nquad), rv(3,3)
      do k = 1, 3_8*LEN2+1_8
        tp(k) = real(k-1, r64)*2.0_r64*pi/real(3_8*LEN2, r64)
      end do
      sxbd2 = 0.0_r64; swbd2 = 0.0_r64; stangbd2 = 0.0_r64; ss = 0.0_r64
      rv = 0.0_r64
      call line3quadr_3dline(sxt, korder, kpols, umatr, nquad, tgl, wgl, Dgl, &
                             3_8*LEN2, tp, nb2, sxbd2, swbd2, stangbd2, ss, rv)
    end block
    block
      real(r64) :: tp(3*LEN3+1), ss(3*LEN3*nquad), rv(3,3)
      do k = 1, 3_8*LEN3+1_8
        tp(k) = real(k-1, r64)*2.0_r64*pi/real(3_8*LEN3, r64)
      end do
      sxbd3 = 0.0_r64; swbd3 = 0.0_r64; stangbd3 = 0.0_r64; ss = 0.0_r64
      rv = 0.0_r64
      call line3quadr_3dline(sxt, korder, kpols, umatr, nquad, tgl, wgl, Dgl, &
                             3_8*LEN3, tp, nb3, sxbd3, swbd3, stangbd3, ss, rv)
    end block

    dl = w_bclag/(-1.0_r64 - tgl)
    dr = w_bclag/( 1.0_r64 - tgl)
    bclagmatlr(:,1) = dl/sum(dl)
    bclagmatlr(:,2) = dr/sum(dr)
    call cpu_time(t1c);  timeinfo(2) = timeinfo(2) + (t1c - t0)

    call cpu_time(t0)
    idxs = 0_8;  ms = 0_8
    call lap3dsdlpmat_levels_r64(m, tx, n, lv, isimd, As, Ad, idxs, ms)
    call cpu_time(t1c);  timeinfo(3) = timeinfo(3) + (t1c - t0)

    if (ms > 0_8) then
      block
        real(r64)    :: txn(3,ms)
        real(r64)    :: w1(nquad,3,ms), w3(nquad,3,ms), w5(nquad,3,ms)
        complex(r64) :: troot(ms,3), xr(ms,3), yr(ms,3), zr(ms,3)
        real(r64)    :: xrr(ms,3), yrr(ms,3), zrr(ms,3)
        integer(8)   :: rfc(ms,3), rfc_ssq(ms), idxnp(ncoeff)
        complex(r64) :: Fx(3*nquad,ncoeff), Fy(3*nquad,ncoeff)
        complex(r64) :: Fz(3*nquad,ncoeff)
        complex(r64) :: Ichi(ms,ncoeff,4), Ialpha(ms,ncoeff,4)
        real(r64)    :: Ias(ms), om_slp(h_dim,ms), om(h_dim,ms,4)
        real(r64)    :: Ocl(h_dim,ms), Oadd(4*h_dim,ms), Atmp(4*h_dim,ms)
        real(r64)    :: MD(4*h_dim,4*h_dim), Msl(h_dim,h_dim)
        real(r64)    :: Fm(h_dim,h_dim), Otmp(h_dim,ms)
        integer(8)   :: ipiv4(4*h_dim), ipiv1(h_dim), info

        do j = 1, ms
          txn(:,j) = alpha*matmul(R, tx(:,idxs(j)) - c)
        end do

        call cpu_time(t0)
        w1 = 0.0_r64;  w3 = 0.0_r64;  w5 = 0.0_r64
        troot = (0.0_r64,0.0_r64);  xr = troot;  yr = troot;  zr = troot
        rfc = 0_8;  rfc_ssq = 0_8
        call build_ssq_weights_r64(ms, txn, RHO_SSQ, nbd, sxbd, sbdnp, nquad, &
                                   tgl, wgl, Legmat, w1, w3, w5, &
                                   troot, xr, yr, zr, rfc, rfc_ssq)
        call cpu_time(t1c);  timeinfo(4) = timeinfo(4) + (t1c - t0)

        call cpu_time(t0)
        xrr = 0.0_r64;  yrr = 0.0_r64;  zrr = 0.0_r64
        xrr = real(xr, r64)
        yrr = real(yr, r64)
        zrr = real(zr, r64)
        Ias = 0.0_r64
        call evaluate_solid_angle_integral_fast_r64(ms, txn, nbd, sbdnp, nquad, &
             sxbd, stangbd, sspbd, &
             LEN1, sxbd1, stangbd1, swbd1, &
             LEN2, sxbd2, stangbd2, swbd2, &
             LEN3, sxbd3, stangbd3, swbd3, &
             qhat, tgl, wgl, Dgl, w_bclag, bclagmatlr, &
             troot, xrr, yrr, zrr, rfc, Ias, &
             rho_in=8.0_r64**(8.0_r64/real(nquad, r64)), &
             sxbd_raw=sxbd_raw, tx_raw=tx(:,idxs(1:ms)), &
             Rfr=R, alpha_fr=alpha, Legmat=Legmat)
        call cpu_time(t1c);  timeinfo(5) = timeinfo(5) + (t1c - t0)

        call cpu_time(t0)
        t1 = 0_8;  t2 = 0_8
        do ij = 0, nterms
          do k = -ij, ij
            t2 = t2 + 1_8
            if (k <= 0_8) then
              t1 = t1 + 1_8;  idxnp(t1) = t2
            end if
          end do
        end do
        Fx = Fxbd(:,idxnp);  Fy = Fybd(:,idxnp);  Fz = Fzbd(:,idxnp)

        call build_tc_chialpha_r64(ms, txn, nbd, sbdnp, nquad, ncoeff, &
                                   sxbd, stangbd, sspbd, w1, w3, &
                                   Fx, Fy, Fz, Ias, Ichi, Ialpha)
        call cpu_time(t1c);  timeinfo(6) = timeinfo(6) + (t1c - t0)

        call cpu_time(t0)
        call qao_omegasdlp_r64(ms, nterms, ncoeff, h_dim, txn, Ichi, Ialpha, &
                               om_slp, om)
        call cpu_time(t1c);  timeinfo(7) = timeinfo(7) + (t1c - t0)

        call cpu_time(t0)
        fac = 1.0_r64/(4.0_r64*pi)/alpha
        Ocl = fac*om_slp
        Oadd(        1:  h_dim, :) =  om(:,:,1)
        Oadd(  h_dim+1:2*h_dim, :) = -om(:,:,2)
        Oadd(2*h_dim+1:3*h_dim, :) = -om(:,:,3)
        Oadd(3*h_dim+1:4*h_dim, :) = -om(:,:,4)
        Oadd = -fac*Oadd

        Msl = Mmatrix(      1:  h_dim,       1:  h_dim)
        Fm = Mmatrix(h_dim+1:2*h_dim, h_dim+1:2*h_dim)
        MD = Mmatrix
        MD(      1:  h_dim,       1:  h_dim) = 0.0_r64
        MD(h_dim+1:2*h_dim, h_dim+1:2*h_dim) = 0.0_r64

        Atmp = Oadd
        MD   = transpose(MD)
        call dgesv(4_8*h_dim, ms, MD, 4_8*h_dim, ipiv4, Atmp, 4_8*h_dim, info)

        Otmp = Ocl + matmul(transpose(Fm), Atmp(1:h_dim,:))
        Msl  = transpose(Msl)
        call dgesv(h_dim, ms, Msl, h_dim, ipiv1, Otmp, h_dim, info)

        do j = 1, ms
          As(idxs(j),:)           = Otmp(:,j)
          Ad(idxs(j),:)           = alpha*Atmp(1:h_dim,j)
          Omega(:,idxs(j))        = Oadd(:,j)
          IalphaAsvestas(idxs(j)) = Ias(j)
        end do
        call cpu_time(t1c);  timeinfo(8) = timeinfo(8) + (t1c - t0)
      end block
    end if

    call clear_patch_levels_r64(lv)

  end subroutine rrq_r64

end module lap3d_close_mod
