
program stellarator_grf
  use iso_c_binding,         only: c_long_long, c_double
  use quatapproximation_mod, only: r64
  use linequaaadrature_mod,  only: gauss_r64
  use koorn_geom_mod,        only: get_vioreanu_nodes, get_vioreanu_wts
  use lap3d_close_mod,       only: rrq_r64
  use lap3d_mod,             only: lap3dsdlpmat_r64
  use patch_refine_mod,      only: PR_ADAPTIVE, PR_SQUARE
  use omp_lib,               only: omp_get_max_threads, omp_get_wtime
  implicit none

  interface
    function stellarator_geo_ntri(mp, np, p, x0) &
             bind(C, name='stellarator_geo_ntri')
      import :: c_long_long, c_double
      integer(c_long_long)             :: stellarator_geo_ntri
      integer(c_long_long), intent(in) :: mp, np, p
      real(c_double),       intent(in) :: x0(*)
    end function stellarator_geo_ntri

    subroutine stellarator_geo_uv2x(mp, np, p, x0, itri, nuv, uv, x) &
               bind(C, name='stellarator_geo_uv2x')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: mp, np, p, itri, nuv
      real(c_double),       intent(in)  :: x0(*), uv(*)
      real(c_double),       intent(out) :: x(*)
    end subroutine stellarator_geo_uv2x

    subroutine stellarator_geo(mp, np, p, x0, D, uvs, wts, &
                               sx, snx, sw, rts, rps) &
               bind(C, name='stellarator_geo')
      import :: c_long_long, c_double
      integer(c_long_long), intent(in)  :: mp, np, p
      real(c_double),       intent(in)  :: x0(*), D(*), uvs(*), wts(*)
      real(c_double),       intent(out) :: sx(*), snx(*), sw(*), rts(*), rps(*)
    end subroutine stellarator_geo

    subroutine lfmm3d_t_cd_p(eps, nsource, source, charge, dipvec, &
                             ntarg, targ, pottarg, ier)
      integer(8) :: nsource, ntarg, ier
      real(8)    :: eps, source(3,*), charge(*), dipvec(3,*)
      real(8)    :: targ(3,*), pottarg(*)
    end subroutine lfmm3d_t_cd_p
  end interface

  real(r64),  parameter :: LAM = 0.5_r64
  integer(8), parameter :: NC = 7_8
  integer(8), parameter :: ORDER_L(NC) = [ 4_8,  6_8,  8_8, 10_8, 12_8, 14_8, 16_8]
  integer(8), parameter :: MP_L(NC)    = [24_8, 24_8, 36_8, 36_8, 48_8, 60_8, 72_8]
  integer(8), parameter :: NP_L(NC)    = [72_8, 72_8,108_8,108_8,144_8,180_8,216_8]
  real(r64),  parameter :: TOL_L(NC)   = [1.0e-04_r64, 1.0e-06_r64, 1.0e-08_r64, &
                                          1.0e-09_r64, 1.0e-11_r64, 1.0e-14_r64, &
                                          1.0e-14_r64]

  integer(8) :: ncases, icase, isimd, ichart, iadap
  integer(8) :: mp, np, nterms, hdim, nquad, orderff, ntri, nsrc
  integer(8) :: k, i, j, ier, nbd, sbdnp, nnz, mtcmax
  real(r64)  :: distff, timeinfo(20), err, fmm_eps
  real(r64)  :: pi, tt, th1, th2, thi1
  real(r64)  :: t_fmm, t_close, t_geo, w0, w1p, w2p, w3p
  integer(8) :: c0, c1, crate
  logical    :: exterior

  real(r64), allocatable :: tgl(:), wgl(:), Dgl(:,:), uvs(:,:), wts(:)
  real(r64), allocatable :: tglq(:), wglq(:), Dglq(:,:), uvbd(:,:), tpan(:)
  real(r64), allocatable :: sx(:,:), snx(:,:), sw(:), rts(:,:), rps(:,:)
  real(r64), allocatable :: ub(:), ubn(:), u(:), charge(:), dipvec(:,:)
  real(r64), allocatable :: csr_val(:)
  integer(8),allocatable :: mtcs(:), offs(:), csr_idx(:)

  character(len=32) :: arg

  ncases = 3_8;  isimd = 0_8;  ichart = 1_8;  iadap = 0_8
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg);  read(arg,*) ncases
  end if
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg);  read(arg,*) isimd
  end if
  if (command_argument_count() >= 3) then
    call get_command_argument(3, arg);  read(arg,*) ichart
  end if
  if (command_argument_count() >= 4) then
    call get_command_argument(4, arg);  read(arg,*) iadap
  end if
  PR_ADAPTIVE = (iadap == 1_8)
  if (iadap == 1_8 .or. iadap == 3_8) PR_SQUARE = .false.

  exterior = .true.
  distff   = 1.4_r64
  pi       = 4.0_r64*atan(1.0_r64)
  write(*,'(a,i0,a)') 'close loop on ', omp_get_max_threads(), &
      ' threads (timeinfo is cpu-summed across threads)'

  do icase = 1, min(ncases, NC)
    nterms  = ORDER_L(icase)
    mp      = MP_L(icase)
    np      = NP_L(icase)
    fmm_eps = TOL_L(icase)
    hdim    = nterms*(nterms+1_8)/2_8
    nquad   = nterms + 2_8
    orderff = nterms + 2_8
    timeinfo = 0.0_r64

    call system_clock(c0, crate)
    allocate(tgl(nterms), wgl(nterms), Dgl(nterms,nterms))
    allocate(uvs(2,hdim), wts(hdim))
    call gauss_r64(nterms, tgl, wgl, Dgl)
    call get_vioreanu_nodes(nterms-1_8, hdim, uvs)
    call get_vioreanu_wts(nterms-1_8, hdim, wts)
    ntri = stellarator_geo_ntri(int(mp,c_long_long), int(np,c_long_long), &
                                int(nterms,c_long_long), tgl)
    nsrc = ntri*hdim
    allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
    call stellarator_geo(int(mp,c_long_long), int(np,c_long_long), &
                         int(nterms,c_long_long), tgl, Dgl, uvs, wts, &
                         sx, snx, sw, rts, rps)
    call system_clock(c1);  t_geo = real(c1-c0, r64)/real(crate, r64)

    allocate(ub(nsrc), ubn(nsrc))
    do i = 1, nsrc
      ub(i)  = f_harm(sx(:,i))
      ubn(i) = dot_product(snx(:,i), gradf_harm(sx(:,i)))
    end do

    allocate(u(nsrc), charge(nsrc), dipvec(3,nsrc))
    charge = sw*ubn
    do i = 1, nsrc
      dipvec(:,i) = -sw(i)*ub(i)*snx(:,i)
    end do
    call system_clock(c0)
    ier = 0_8
    call lfmm3d_t_cd_p(fmm_eps, nsrc, sx, charge, dipvec, nsrc, sx, u, ier)
    call system_clock(c1);  t_fmm = real(c1-c0, r64)/real(crate, r64)
    if (ier /= 0_8) write(*,'(a,i0)') 'WARNING: lfmm3d ier = ', ier

    sbdnp = 3_8;  nbd = sbdnp*nquad
    allocate(tglq(nquad), wglq(nquad), Dglq(nquad,nquad))
    allocate(tpan(sbdnp+1), uvbd(2,nbd+3))
    call gauss_r64(nquad, tglq, wglq, Dglq)
    do k = 1, sbdnp+1_8
      tpan(k) = real(k-1, r64)*2.0_r64*pi/real(sbdnp, r64)
    end do
    th1 = 2.0_r64*pi/3.0_r64;  th2 = 4.0_r64*pi/3.0_r64;  thi1 = 1.0_r64/th1
    do k = 1, sbdnp
      do i = 1, nquad
        j  = (k-1_8)*nquad + i
        tt = tpan(k) + 0.5_r64*(1.0_r64 + tglq(i))*(tpan(k+1) - tpan(k))
        uvbd(:,j) = tpar2uv(tt, th1, th2, thi1)
      end do
    end do
    uvbd(:,nbd+1) = [0.0_r64, 0.0_r64]
    uvbd(:,nbd+2) = [1.0_r64, 0.0_r64]
    uvbd(:,nbd+3) = [0.0_r64, 1.0_r64]

    allocate(mtcs(ntri), offs(ntri+1_8))
    call system_clock(c0)
    w0 = omp_get_wtime()

    !$omp parallel do schedule(static)
    do k = 1, ntri
      block
        real(r64) :: qr, qp(3), swk(hdim), sxk(3,hdim)
        integer(8) :: i2, j2, cnt
        do i2 = 1, hdim
          sxk(:,i2) = sx(:, (k-1_8)*hdim + i2)
          swk(i2)   = sw((k-1_8)*hdim + i2)
        end do
        qr    = 1.75_r64*sqrt(sum(swk))
        qp(1) = sum(sxk(1,:))/real(hdim, r64)
        qp(2) = sum(sxk(2,:))/real(hdim, r64)
        qp(3) = sum(sxk(3,:))/real(hdim, r64)
        cnt = 0_8
        do j2 = 1, nsrc
          if ((sx(1,j2)-qp(1))**2 + (sx(2,j2)-qp(2))**2 + &
              (sx(3,j2)-qp(3))**2 < qr**2) cnt = cnt + 1_8
        end do
        mtcs(k) = cnt
      end block
    end do
    !$omp end parallel do
    offs(1) = 0_8
    do k = 1, ntri
      offs(k+1_8) = offs(k) + mtcs(k)
    end do
    nnz = offs(ntri+1_8)
    allocate(csr_idx(nnz), csr_val(nnz))
    w1p = omp_get_wtime()

    block
      real(r64) :: sxk(3,hdim), snxk(3,hdim), swk(hdim)
      real(r64) :: rtsk(3,hdim), rpsk(3,hdim), xbuf(3,nbd+3)
      real(r64) :: tw(3,1), S1w(1,hdim), K1w(1,hdim), O1w(4*hdim,1), I1w(1)
      real(r64) :: tinfo0(20)
      integer(8) :: i2
      do i2 = 1, hdim
        sxk(:,i2)  = sx(:,i2);   snxk(:,i2) = snx(:,i2);  swk(i2) = sw(i2)
        rtsk(:,i2) = rts(:,i2);  rpsk(:,i2) = rps(:,i2)
      end do
      tw(:,1) = sx(:,1) + 0.01_r64*snx(:,1)
      xbuf = 0.0_r64
      if (ichart == 1_8) &
        call stellarator_geo_uv2x(int(mp,c_long_long), int(np,c_long_long), &
               int(nterms,c_long_long), tgl, int(1,c_long_long), &
               int(nbd+3,c_long_long), uvbd, xbuf)
      S1w = 0.0_r64;  K1w = 0.0_r64;  O1w = 0.0_r64;  I1w = 0.0_r64
      tinfo0 = 0.0_r64
      call rrq_r64(1_8, tw, hdim, sxk, snxk, swk, rtsk, rpsk, &
                   nterms, nquad, orderff, distff, exterior, isimd, &
                   ichart, xbuf(:,1:nbd), xbuf(:,nbd+1:nbd+3), &
                   S1w, K1w, O1w, I1w, tinfo0)
    end block

    mtcmax = maxval(mtcs)
    !$omp parallel default(shared) reduction(+:timeinfo)
    block
      real(r64), allocatable :: wtcx(:), wAs(:), wAd(:)
      real(r64), allocatable :: wS(:), wK(:), wOm(:), wIa(:)
      allocate(wtcx(3*mtcmax), wAs(mtcmax*hdim), wAd(mtcmax*hdim))
      allocate(wS(mtcmax*hdim), wK(mtcmax*hdim))
      allocate(wOm(4*hdim*mtcmax), wIa(mtcmax))

    !$omp do schedule(dynamic,4)
    do k = 1, ntri
      if (mtcs(k) /= 0_8) then
        block
          real(r64) :: qr, qp(3), xbuf(3,nbd+3)
          real(r64) :: sxk(3,hdim), snxk(3,hdim), swk(hdim)
          real(r64) :: rtsk(3,hdim), rpsk(3,hdim)
          integer(8) :: i2, j2, mtc, o0, i0
          i0 = (k-1_8)*hdim
          do i2 = 1, hdim
            sxk(:,i2)  = sx(:,i0+i2);   snxk(:,i2) = snx(:,i0+i2)
            swk(i2)    = sw(i0+i2)
            rtsk(:,i2) = rts(:,i0+i2);  rpsk(:,i2) = rps(:,i0+i2)
          end do
          qr    = 1.75_r64*sqrt(sum(swk))
          qp(1) = sum(sxk(1,:))/real(hdim, r64)
          qp(2) = sum(sxk(2,:))/real(hdim, r64)
          qp(3) = sum(sxk(3,:))/real(hdim, r64)
          mtc = mtcs(k);  o0 = offs(k)
          j2 = 0_8
          do i2 = 1, nsrc
            if ((sx(1,i2)-qp(1))**2 + (sx(2,i2)-qp(2))**2 + &
                (sx(3,i2)-qp(3))**2 < qr**2) then
              j2 = j2 + 1_8
              if (j2 > mtc) exit
              csr_idx(o0+j2) = i2
            end if
          end do
          if (j2 /= mtc) then
            write(*,'(a,i0,a,i0,a,i0)') 'COUNT MISMATCH panel ', k, &
                ': phase0 ', mtc, '  phase1 ', j2
            error stop 'csr count mismatch'
          end if

          do j2 = 1, mtc
            wtcx(3*j2-2_8:3*j2) = sx(:, csr_idx(o0+j2))
          end do
          wS(1:mtc*hdim)     = 0.0_r64
          wK(1:mtc*hdim)     = 0.0_r64
          wOm(1:4*hdim*mtc)  = 0.0_r64
          wIa(1:mtc)         = 0.0_r64
          xbuf = 0.0_r64
          if (ichart == 1_8) &
            call stellarator_geo_uv2x(int(mp,c_long_long), int(np,c_long_long), &
                   int(nterms,c_long_long), tgl, int(k,c_long_long), &
                   int(nbd+3,c_long_long), uvbd, xbuf)
          call rrq_r64(mtc, wtcx, hdim, sxk, snxk, swk, rtsk, rpsk, &
                       nterms, nquad, orderff, distff, exterior, isimd, &
                       ichart, xbuf(:,1:nbd), xbuf(:,nbd+1:nbd+3), &
                       wS, wK, wOm, wIa, timeinfo)

          call lap3dsdlpmat_r64(mtc, wtcx, hdim, sxk, snxk, swk, wAs, wAd)

          do j2 = 1, mtc
            csr_val(o0+j2) = 0.0_r64
          end do
          do i2 = 1, hdim
            do j2 = 1, mtc
              csr_val(o0+j2) = csr_val(o0+j2) &
                + (wS((i2-1_8)*mtc+j2) - wAs((i2-1_8)*mtc+j2))*ubn(i0+i2) &
                - (wK((i2-1_8)*mtc+j2) - wAd((i2-1_8)*mtc+j2))*ub(i0+i2)
            end do
          end do
        end block
      end if
    end do
    !$omp end do
      deallocate(wtcx, wAs, wAd, wS, wK, wOm, wIa)
    end block
    !$omp end parallel
    w2p = omp_get_wtime()

    do k = 1, ntri
      do j = offs(k)+1_8, offs(k+1_8)
        u(csr_idx(j)) = u(csr_idx(j)) + csr_val(j)
      end do
    end do
    w3p = omp_get_wtime()
    call system_clock(c1);  t_close = real(c1-c0, r64)/real(crate, r64)
    write(*,'(a,f8.3,a,f8.3,a,f8.3,a)') '  phases: ballcount ', w1p-w0, &
        ' s | rrq+fill ', w2p-w1p, ' s | scatter ', w3p-w2p, ' s'

    err = maxval(abs(u))/maxval(abs(ubn))
    write(*,'(a,i0,a,i2,a,i3,a,i3,a,i0,a,es10.3,a,f8.2,a,f8.2,a)') &
        'case ', icase, ': order=', nterms, '  mp=', mp, '  np=', np, &
        '  N=', nsrc, '  GRF max rel err = ', err, &
        '  (fmm ', t_fmm, ' s, close ', t_close, ' s)'
    write(*,'(a,8es10.2)') '  timeinfo ', timeinfo(1:8)

    deallocate(tgl, wgl, Dgl, uvs, wts, tglq, wglq, Dglq, tpan, uvbd)
    deallocate(sx, snx, sw, rts, rps, ub, ubn, u, charge, dipvec)
    deallocate(mtcs, offs, csr_idx, csr_val)
  end do

contains

  pure function tpar2uv(t, thold1, thold2, tholdinv1) result(uv)
    real(r64), intent(in) :: t, thold1, thold2, tholdinv1
    real(r64) :: uv(2)
    uv = 0.0_r64
    if (t < thold1) then
      uv(2) = 1.0_r64 - tholdinv1*t
    else if (t < thold2) then
      uv(1) = tholdinv1*(t - thold1)
    else
      uv(1) = 1.0_r64 - tholdinv1*(t - thold2)
      uv(2) = tholdinv1*(t - thold2)
    end if
  end function tpar2uv

  pure function f_harm(x) result(v)
    real(r64), intent(in) :: x(3)
    real(r64) :: v
    v = exp(LAM*x(1))*cos(LAM*x(2))
  end function f_harm

  pure function gradf_harm(x) result(g)
    real(r64), intent(in) :: x(3)
    real(r64) :: g(3)
    g(1) =  LAM*exp(LAM*x(1))*cos(LAM*x(2))
    g(2) = -LAM*exp(LAM*x(1))*sin(LAM*x(2))
    g(3) =  0.0_r64
  end function gradf_harm

end program stellarator_grf
