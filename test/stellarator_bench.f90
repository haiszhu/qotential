
program stellarator_bench
  use quatapproximation_mod, only: r64
  use lap3d_close_mod, only: rrq_r64
  implicit none

  character(len=200) :: DATDIR

  integer(8), parameter :: NC = 6_8
  integer(8), parameter :: ORDER_L(NC) = [ 4_8,  6_8,  8_8, 10_8, 12_8, 14_8]
  integer(8), parameter :: MP_L(NC)    = [24_8, 24_8, 36_8, 36_8, 48_8, 60_8]
  integer(8), parameter :: NP_L(NC)    = [72_8, 72_8,108_8,108_8,144_8,180_8]

  integer(8) :: nterms, mp, np, npan, nsrc, m, k, j, i, hdim, nquad, orderff
  integer(8) :: ncases, icase, itarg, isimd
  real(r64)  :: mpr, npr, ntr, hdr, nqr, npr2, ofr
  real(r64)  :: distff, timeinfo(20), qradii, qpoint(3), tsrc, ttgt
  logical    :: exterior

  real(r64),  allocatable :: sx(:,:), snx(:,:), sw(:), rts(:,:), rps(:,:)
  real(r64),  allocatable :: sxk(:,:), snxk(:,:), swk(:), rtsk(:,:), rpsk(:,:)
  real(r64),  allocatable :: tcx(:,:), S_ij(:,:), K_ij(:,:)
  real(r64),  allocatable :: Omegas(:,:), IalphaAsvestas(:)
  logical,    allocatable :: near(:)
  integer(8), allocatable :: idxk(:)

  character(len=200) :: fname
  character(len=32)  :: arg

  ncases = 1_8
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg);  read(arg,*) ncases
  end if
  isimd = 0_8
  if (command_argument_count() >= 2) then
    call get_command_argument(2, arg);  read(arg,*) isimd
  end if
  DATDIR = './'
  if (command_argument_count() >= 3) call get_command_argument(3, DATDIR)
  exterior = .true.          ! rrq's iside = 30

  do icase = 1, min(ncases, NC)
    nterms = ORDER_L(icase)
    mp     = MP_L(icase)
    np     = NP_L(icase)

    write(fname,'(a,a,i0,a,i0,a,i0,a)') trim(DATDIR), &
        'stellarator_source_mp', mp, '_np', np, '_order', nterms, '.dat'
    open(1, file=trim(fname), status='old')
    read(1,*) mpr, npr, ntr, hdr, nqr, npr2, ofr
    hdim    = int(hdr,  8)
    nquad   = int(nqr,  8)
    npan    = int(npr2, 8)
    orderff = int(ofr,  8)
    nsrc    = npan*nterms*(nterms+1_8)/2_8
    allocate(sx(3,nsrc), snx(3,nsrc), sw(nsrc), rts(3,nsrc), rps(3,nsrc))
    do i = 1, nsrc
      read(1,*) sx(1:3,i), snx(1:3,i), sw(i), rts(1:3,i), rps(1:3,i)
    end do
    close(1)

    distff   = 1.4_r64
    timeinfo = 0.0_r64
    allocate(idxk(hdim), sxk(3,hdim), snxk(3,hdim), swk(hdim))
    allocate(rtsk(3,hdim), rpsk(3,hdim), near(nsrc))

    do k = 1, npan
      do i = 1, hdim
        idxk(i) = (k-1_8)*hdim + i
      end do
      sxk  = sx(:,idxk);   swk  = sw(idxk);   snxk = snx(:,idxk)
      rtsk = rts(:,idxk);  rpsk = rps(:,idxk)

      qradii    = 1.75_r64*sqrt(sum(swk))
      qpoint(1) = sum(sxk(1,:))/real(hdim, r64)
      qpoint(2) = sum(sxk(2,:))/real(hdim, r64)
      qpoint(3) = sum(sxk(3,:))/real(hdim, r64)
      do j = 1, nsrc
        near(j) = sum((sx(:,j) - qpoint)**2) < qradii**2
      end do
      m = count(near)

      allocate(tcx(3,m), S_ij(m,hdim), K_ij(m,hdim))
      allocate(Omegas(4*hdim,m), IalphaAsvestas(m))
      itarg = 0_8
      do j = 1, nsrc
        if (near(j)) then
          itarg = itarg + 1_8;  tcx(:,itarg) = sx(:,j)
        end if
      end do

      call rrq_r64(m, tcx, hdim, sxk, snxk, swk, rtsk, rpsk, &
                   nterms, nquad, orderff, distff, exterior, isimd, &
                   S_ij, K_ij, Omegas, IalphaAsvestas, timeinfo)

      deallocate(tcx, S_ij, K_ij, Omegas, IalphaAsvestas)
    end do

    tsrc = sum(timeinfo(1:2))
    ttgt = sum(timeinfo(3:8))
    write(*,'(a)') '======================'
    write(*,'(a,i4,a,i4,a,i4,a,i4)') 'Case: nterms =', nterms, &
        ',      mp =', mp, ',      np =', np, ',   isimd =', isimd
    write(*,'(a,i0)')      ' nsrc        ', nsrc
    write(*,'(a,es12.4)')  ' source time ', tsrc
    write(*,'(a,es12.4)')  ' source pps  ', real(nsrc,r64)/tsrc
    write(*,'(a,es12.4)')  ' target time ', ttgt
    write(*,'(a,es12.4)')  ' target pps  ', real(nsrc,r64)/ttgt
    write(*,'(a,8es10.2)') ' timeinfo    ', timeinfo(1:8)
    write(*,'(a)') ' '

    deallocate(sx, snx, sw, rts, rps)
    deallocate(idxk, sxk, snxk, swk, rtsk, rpsk, near)
  end do

end program stellarator_bench
