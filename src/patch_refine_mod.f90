module patch_refine_mod

  use quatapproximation_mod, only: r64
  use koorn_geom_mod, only: koorn_vals2coefs_coefs2vals, koorn_pols_batch_r64, &
                            get_vioreanu_nodes, get_vioreanu_wts
  use lap3d_mod, only: lap3dsdlpmat_r64
  use lap3d_simd_mod, only: csimd128lap3dsdlpmat_r64, &
                            csimd256lap3dsdlpmat_r64, &
                            csimd512lap3dsdlpmat_r64

  implicit none
  private
  public :: patch_level_t, patch_levels_t, build_patch_levels_r64, &
            lap3dsdlpmat_levels_r64, patch_ref_t, patch_ref_level_t, REF, &
            PR_ADAPTIVE, PR_SQUARE, PR_QGL

  integer(8), parameter :: CC_SUBIDX(39) = [ &
     1_8,  2_8,  3_8,  4_8,  7_8,  8_8,  9_8, 10_8, 11_8, 13_8, &
    14_8, 16_8, 17_8, 27_8, 32_8, 33_8, 34_8, 35_8, 37_8, 40_8, &
    41_8, 42_8, 43_8, 44_8, 46_8, 47_8, 48_8, 49_8, 50_8, 52_8, &
    53_8, 55_8, 58_8, 59_8, 60_8, 61_8, 62_8, 63_8, 64_8 ]

  real(r64), parameter :: CC_EXTRA_SCALE  = 0.625_r64
  real(r64), parameter :: CC_EXTRA_OFFSET = 0.4375_r64

  type :: patch_level_t
    integer(8) :: n    = 0_8            ! quadrature nodes on this level
    integer(8) :: nsub = 0_8            ! sub-patches on this level
    real(r64), allocatable :: sx(:,:)        ! (3,n)  surface points
    real(r64), allocatable :: snx(:,:)       ! (3,n)  unit normals
    real(r64), allocatable :: sw(:)          ! (n)    area weights
    real(r64), allocatable :: qpoint(:,:)    ! (3,nsub) weighted centroid
    real(r64), allocatable :: qradii2(:)     ! (nsub)   squared near radius
  end type patch_level_t

  type :: patch_ref_level_t
    integer(8) :: n = 0_8, nsub = 0_8
    real(r64), allocatable :: uv(:,:), wt(:), interpmat(:,:)
    integer(8), allocatable :: suboff(:), subtyp(:)
  end type patch_ref_level_t

  type :: patch_ref_t
    integer(8) :: korder = -1_8, npols = -1_8, order = -1_8, nlevel = -1_8
    logical    :: sqmode = .false.
    integer(8) :: qgl = -1_8
    type(patch_ref_level_t), allocatable :: lev(:)
  end type patch_ref_t

  type(patch_ref_t), save :: REF

  logical, save :: PR_ADAPTIVE = .false.

  logical,    save :: PR_SQUARE = .true.
  integer(8), save :: PR_QGL = 0_8
  real(r64), allocatable, save :: WS_Ab(:,:), WS_Bb(:,:), WS_Cb(:,:)
  !$omp threadprivate(WS_Ab, WS_Bb, WS_Cb)

  type :: cut_cache_t
    integer(8) :: mask = -1_8, ng = 0_8
    real(r64), allocatable :: gI(:,:)
  end type cut_cache_t
  integer(8), parameter :: NCUTC = 128_8
  type(cut_cache_t), save :: CUTC(NCUTC)
  integer(8), save :: ncutn = 0_8, cutc_key = -1_8
  !$omp threadprivate(CUTC, ncutn, cutc_key)

  type :: patch_levels_t
    integer(8) :: nlevel = 0_8
    type(patch_level_t), allocatable :: lev(:)   ! 1=ff, 2=f, 3=c, 4=cc
  end type patch_levels_t

contains

  subroutine build_patch_ref_r64(korder, npols, order, nlevel)
    integer(8), intent(in) :: korder, npols, order, nlevel

    real(r64)  :: umatr(npols,npols), vmatr(npols,npols)
    integer(8) :: ilev, n1, nprev, n, nsub, ncc, i_off

    if (REF%korder == korder .and. REF%npols  == npols .and. &
        REF%order  == order  .and. REF%nlevel == nlevel) return

    block
      integer(8) :: qeff
      qeff = -1_8
      if (PR_SQUARE) then
        qeff = korder + 1_8
        if (PR_QGL > 0_8) qeff = PR_QGL
      end if
      if (PR_SQUARE .and. PR_ADAPTIVE) &
        error stop 'PR_SQUARE and PR_ADAPTIVE cannot be combined'
    if (.not. (REF%korder == korder .and. REF%npols  == npols .and. &
               REF%order  == order  .and. REF%nlevel == nlevel .and. &
               (REF%sqmode .eqv. PR_SQUARE) .and. REF%qgl == qeff)) then

    if (allocated(REF%lev)) deallocate(REF%lev)
    allocate(REF%lev(nlevel))

    umatr = 0.0_r64;  vmatr = 0.0_r64
    call koorn_vals2coefs_coefs2vals(korder, npols, umatr, vmatr)

    n1  = order*(order+1_8)/2_8        ! nodes per sub-patch, VR order-1
    ncc = size(CC_SUBIDX, kind=8) + 4_8

    do ilev = 1, nlevel

      if (ilev == 1_8) then
        n = n1;  nsub = 1_8
        allocate(REF%lev(1)%uv(2,n), REF%lev(1)%wt(n))
        call get_vioreanu_nodes(order-1_8, n, REF%lev(1)%uv)
        call get_vioreanu_wts  (order-1_8, n, REF%lev(1)%wt)

      else if (ilev == nlevel .and. nlevel == 4_8) then
        n = ncc*n1;  nsub = ncc
        allocate(REF%lev(ilev)%uv(2,n), REF%lev(ilev)%wt(n))
        block
          real(r64), allocatable :: uvt(:,:), wtt(:)
          real(r64)  :: uv0(2)
          integer(8) :: np, q, b, i, isrc, idst
          np = REF%lev(ilev-1)%n
          allocate(uvt(2,4*np), wtt(4*np))
          do i = 1, np
            uv0 = 2.0_r64*REF%lev(ilev-1)%uv(:,i) - 1.0_r64
            uvt(:,      i) =  uv0
            uvt(:,   np+i) = -uv0
            uvt(1, 2*np+i) =  uv0(1)
            uvt(2, 2*np+i) =  uv0(2) + 2.0_r64
            uvt(1, 3*np+i) =  uv0(1) + 2.0_r64
            uvt(2, 3*np+i) =  uv0(2)
            do b = 0, 3
              wtt(b*np+i) = 0.25_r64*REF%lev(ilev-1)%wt(i)
            end do
          end do
          uvt = ((uvt - 1.0_r64)*0.5_r64 + 1.0_r64)*0.5_r64
          do q = 1, size(CC_SUBIDX, kind=8)
            do i = 1, n1
              isrc = (CC_SUBIDX(q)-1_8)*n1 + i
              idst = (q-1_8)*n1 + i
              REF%lev(ilev)%uv(:,idst) = uvt(:,isrc)
              REF%lev(ilev)%wt(idst)   = wtt(isrc)
            end do
          end do
          do q = 1, 4
            do i = 1, n1
              isrc = (q-1_8)*n1 + i
              idst = (size(CC_SUBIDX, kind=8) + q - 1_8)*n1 + i
              REF%lev(ilev)%uv(:,idst) = &
                   (REF%lev(2)%uv(:,isrc) - 0.5_r64)*CC_EXTRA_SCALE + CC_EXTRA_OFFSET
              REF%lev(ilev)%wt(idst)   = CC_EXTRA_SCALE**2 * REF%lev(2)%wt(isrc)
            end do
          end do
          deallocate(uvt, wtt)
        end block

      else
        nprev = REF%lev(ilev-1)%n
        n = 4_8*nprev;  nsub = 4_8*REF%lev(ilev-1)%nsub
        allocate(REF%lev(ilev)%uv(2,n), REF%lev(ilev)%wt(n))
        block
          real(r64)  :: uv0(2)
          integer(8) :: i, b
          do i = 1, nprev
            uv0 = 2.0_r64*REF%lev(ilev-1)%uv(:,i) - 1.0_r64
            REF%lev(ilev)%uv(:,         i) =  uv0
            REF%lev(ilev)%uv(:,  nprev+i) = -uv0
            REF%lev(ilev)%uv(1,2*nprev+i) =  uv0(1)
            REF%lev(ilev)%uv(2,2*nprev+i) =  uv0(2) + 2.0_r64
            REF%lev(ilev)%uv(1,3*nprev+i) =  uv0(1) + 2.0_r64
            REF%lev(ilev)%uv(2,3*nprev+i) =  uv0(2)
            do b = 0, 3
              REF%lev(ilev)%wt(b*nprev+i) = 0.25_r64*REF%lev(ilev-1)%wt(i)
            end do
          end do
          REF%lev(ilev)%uv = ((REF%lev(ilev)%uv - 1.0_r64)*0.5_r64 + 1.0_r64)*0.5_r64
        end block
      end if

      REF%lev(ilev)%n    = n
      REF%lev(ilev)%nsub = nsub

    end do

    do ilev = 1, nlevel
      n    = REF%lev(ilev)%n
      nsub = REF%lev(ilev)%nsub
      if (PR_SQUARE .and. ilev >= 2_8) then
        call ref_square_regroup(ilev, n1, korder, &
             (ilev == nlevel .and. nlevel == 4_8), n, nsub)
        REF%lev(ilev)%n    = n
        REF%lev(ilev)%nsub = nsub
      else
        allocate(REF%lev(ilev)%suboff(nsub+1_8), REF%lev(ilev)%subtyp(nsub))
        REF%lev(ilev)%subtyp = 0_8
        do i_off = 0, nsub
          REF%lev(ilev)%suboff(i_off+1_8) = i_off*n1
        end do
      end if

      allocate(REF%lev(ilev)%interpmat(n,npols))
      block
        real(r64), allocatable :: pols(:,:)
        allocate(pols(n,npols))
        call koorn_pols_batch_r64(n, REF%lev(ilev)%uv, korder, npols, pols)
        call dgemm('N', 'N', n, npols, npols, 1.0_r64, pols, n, &
                   umatr, npols, 0.0_r64, REF%lev(ilev)%interpmat, n)
        deallocate(pols)
      end block

    end do

    REF%korder = korder;  REF%npols  = npols
    REF%order  = order;   REF%nlevel = nlevel
    REF%sqmode = PR_SQUARE
    REF%qgl    = qeff

    end if
    end block

  end subroutine build_patch_ref_r64

  subroutine ref_square_regroup(ilev, n1, korder, is_cc, n, nsub)
    use linequaaadrature_mod, only: gauss_r64
    integer(8), intent(in)    :: ilev, n1, korder
    logical,    intent(in)    :: is_cc
    integer(8), intent(inout) :: n, nsub

    integer(8) :: nsub_old, npairable, ncell, qg, i, j, k, q, ic, jc, row, i0
    integer(8) :: npair, nsingle
    real(r64)  :: side, fac, cu, cv
    real(r64),  allocatable :: uv_old(:,:), wt_old(:)
    real(r64),  allocatable :: tq(:), wq(:), Dq(:,:), u01(:), w01(:)
    integer(8), allocatable :: cellcnt(:,:), singles(:)

    nsub_old  = nsub
    npairable = nsub_old
    if (is_cc) then
      npairable = nsub_old - 4_8
      ncell     = 8_8
    else
      ncell     = 2_8**(ilev-1_8)
    end if
    side = 1.0_r64/real(ncell, r64)
    qg   = korder + 1_8
    if (qg < 8_8) qg = korder + 3_8
    if (PR_QGL > 0_8) qg = PR_QGL
    fac  = 2.0_r64*sum(REF%lev(1)%wt)

    call move_alloc(REF%lev(ilev)%uv, uv_old)
    call move_alloc(REF%lev(ilev)%wt, wt_old)

    allocate(cellcnt(0:ncell-1_8, 0:ncell-1_8), singles(nsub_old))
    cellcnt = 0_8
    do q = 1, npairable
      i0 = (q-1_8)*n1
      cu = sum(uv_old(1, i0+1_8:i0+n1))/real(n1, r64)
      cv = sum(uv_old(2, i0+1_8:i0+n1))/real(n1, r64)
      ic = min(max(int(cu/side, 8), 0_8), ncell-1_8)
      jc = min(max(int(cv/side, 8), 0_8), ncell-1_8)
      cellcnt(ic,jc) = cellcnt(ic,jc) + 1_8
      if (cellcnt(ic,jc) > 2_8) error stop 'square regroup: >2 in a cell'
    end do

    npair = 0_8;  nsingle = 0_8
    do q = 1, npairable
      i0 = (q-1_8)*n1
      cu = sum(uv_old(1, i0+1_8:i0+n1))/real(n1, r64)
      cv = sum(uv_old(2, i0+1_8:i0+n1))/real(n1, r64)
      ic = min(max(int(cu/side, 8), 0_8), ncell-1_8)
      jc = min(max(int(cv/side, 8), 0_8), ncell-1_8)
      if (cellcnt(ic,jc) == 1_8) then
        nsingle = nsingle + 1_8
        singles(nsingle) = q
      end if
    end do
    do q = npairable+1_8, nsub_old
      nsingle = nsingle + 1_8
      singles(nsingle) = q
    end do
    do jc = 0, ncell-1_8
      do ic = 0, ncell-1_8
        if (cellcnt(ic,jc) == 2_8) npair = npair + 1_8
      end do
    end do

    n    = npair*qg*qg + nsingle*n1
    nsub = npair + nsingle
    allocate(REF%lev(ilev)%uv(2,n), REF%lev(ilev)%wt(n))
    allocate(REF%lev(ilev)%suboff(nsub+1_8), REF%lev(ilev)%subtyp(nsub))
    REF%lev(ilev)%subtyp(1:npair)      = 1_8
    REF%lev(ilev)%subtyp(npair+1_8:)   = 0_8

    allocate(tq(qg), wq(qg), Dq(qg,qg), u01(qg), w01(qg))
    call gauss_r64(qg, tq, wq, Dq)
    u01 = 0.5_r64*(tq + 1.0_r64)
    w01 = 0.5_r64*wq

    row = 0_8;  k = 0_8
    REF%lev(ilev)%suboff(1) = 0_8
    do jc = 0, ncell-1_8
      do ic = 0, ncell-1_8
        if (cellcnt(ic,jc) /= 2_8) cycle
        do j = 1, qg
          do i = 1, qg
            row = row + 1_8
            REF%lev(ilev)%uv(1,row) = (real(ic,r64) + u01(i))*side
            REF%lev(ilev)%uv(2,row) = (real(jc,r64) + u01(j))*side
            REF%lev(ilev)%wt(row)   = w01(i)*w01(j)*side*side*fac
          end do
        end do
        k = k + 1_8
        REF%lev(ilev)%suboff(k+1_8) = row
      end do
    end do
    do i = 1, nsingle
      q  = singles(i)
      i0 = (q-1_8)*n1
      do j = 1, n1
        row = row + 1_8
        REF%lev(ilev)%uv(:,row) = uv_old(:, i0+j)
        REF%lev(ilev)%wt(row)   = wt_old(i0+j)
      end do
      k = k + 1_8
      REF%lev(ilev)%suboff(k+1_8) = row
    end do

    deallocate(uv_old, wt_old, cellcnt, singles, tq, wq, Dq, u01, w01)
  end subroutine ref_square_regroup

  subroutine build_patch_levels_r64(korder, npols, sx, rts, rps, order, &
                                    dist, nlevel, lv)
    integer(8),           intent(in)  :: korder, npols, order, nlevel
    real(r64),            intent(in)  :: sx(3,npols), rts(3,npols), rps(3,npols)
    real(r64),            intent(in)  :: dist
    type(patch_levels_t), intent(out) :: lv

    real(r64)  :: srcT(npols,9), dist2
    integer(8) :: ilev, n1, n, nsub

    call build_patch_ref_r64(korder, npols, order, nlevel)

    dist2 = dist*dist

    srcT(:,1:3) = transpose(sx)
    srcT(:,4:6) = transpose(rts)
    srcT(:,7:9) = transpose(rps)

    n1 = order*(order+1_8)/2_8
    lv%nlevel = nlevel
    allocate(lv%lev(nlevel))

    do ilev = 1, nlevel
      n    = REF%lev(ilev)%n
      nsub = REF%lev(ilev)%nsub
      lv%lev(ilev)%n    = n
      lv%lev(ilev)%nsub = nsub

      allocate(lv%lev(ilev)%sx(3,n), lv%lev(ilev)%snx(3,n), lv%lev(ilev)%sw(n))
      allocate(lv%lev(ilev)%qpoint(3,nsub), lv%lev(ilev)%qradii2(nsub))

      block
        real(r64), allocatable :: geo(:,:)
        real(r64)  :: sp, swsum, qp(3), d2, rmax2
        integer(8) :: i, k, i0, nk

        allocate(geo(n,9))
        call dgemm('N', 'N', n, 9_8, npols, 1.0_r64, REF%lev(ilev)%interpmat, n, &
                   srcT, npols, 0.0_r64, geo, n)

        lv%lev(ilev)%sx = transpose(geo(:,1:3))
        do i = 1, n
          lv%lev(ilev)%snx(1,i) = geo(i,8)*geo(i,6) - geo(i,9)*geo(i,5)
          lv%lev(ilev)%snx(2,i) = geo(i,9)*geo(i,4) - geo(i,7)*geo(i,6)
          lv%lev(ilev)%snx(3,i) = geo(i,7)*geo(i,5) - geo(i,8)*geo(i,4)
          sp = sqrt(dot_product(lv%lev(ilev)%snx(:,i), lv%lev(ilev)%snx(:,i)))
          lv%lev(ilev)%snx(:,i) = lv%lev(ilev)%snx(:,i)/sp
          lv%lev(ilev)%sw(i)    = sp*REF%lev(ilev)%wt(i)
        end do

        do k = 1, nsub
          i0    = REF%lev(ilev)%suboff(k)
          nk    = REF%lev(ilev)%suboff(k+1_8) - i0
          swsum = sum(lv%lev(ilev)%sw(i0+1:i0+nk))
          qp(1) = dot_product(lv%lev(ilev)%sx(1,i0+1:i0+nk), &
                              lv%lev(ilev)%sw(i0+1:i0+nk))/swsum
          qp(2) = dot_product(lv%lev(ilev)%sx(2,i0+1:i0+nk), &
                              lv%lev(ilev)%sw(i0+1:i0+nk))/swsum
          qp(3) = dot_product(lv%lev(ilev)%sx(3,i0+1:i0+nk), &
                              lv%lev(ilev)%sw(i0+1:i0+nk))/swsum
          rmax2 = 0.0_r64
          do i = i0+1, i0+nk
            d2 = sum((lv%lev(ilev)%sx(:,i) - qp)**2)
            rmax2 = max(rmax2, d2)
          end do
          lv%lev(ilev)%qpoint(:,k) = qp
          if (REF%lev(ilev)%subtyp(k) == 1_8) then
            lv%lev(ilev)%qradii2(k) = dist2*max(rmax2, 0.5_r64*swsum)
          else
            lv%lev(ilev)%qradii2(k) = dist2*max(rmax2, swsum)
          end if
        end do

        deallocate(geo)
      end block
    end do

  end subroutine build_patch_levels_r64

  subroutine lap3dsdlpmat_levels_r64(m, tx, npols, lv, isimd, As, Ad, idxs, ms)
    integer(8),           intent(in)    :: m, npols, isimd
    real(r64),            intent(in)    :: tx(3,m)
    type(patch_levels_t), intent(in)    :: lv
    real(r64),            intent(inout) :: As(m,npols), Ad(m,npols)
    integer(8),           intent(inout) :: idxs(m), ms

    integer(8) :: ilev, j, k, jb, nrem, nbat
    integer(8) :: irem(m), ibat(m), inext(m)
    logical    :: outside

    nrem = m
    do j = 1, m
      irem(j) = j
    end do

    if (PR_ADAPTIVE .and. lv%nlevel >= 3_8) then
      call adaptive_pass(m, tx, npols, lv, isimd, As, Ad, irem, nrem)
      ilev = lv%nlevel
      call one_level(ilev, m, tx, npols, lv, isimd, As, Ad, irem, nrem)
      ms = nrem
      idxs(1:ms) = irem(1:ms)
      return
    end if

    do ilev = 1, lv%nlevel
      if (nrem == 0_8) exit

      nbat = 0_8
      jb   = 0_8
      do j = 1, nrem
        outside = .true.
        do k = 1, lv%lev(ilev)%nsub
          if (sum((lv%lev(ilev)%qpoint(:,k) - tx(:,irem(j)))**2) &
              < lv%lev(ilev)%qradii2(k)) then
            outside = .false.
            exit
          end if
        end do
        if (outside) then
          nbat = nbat + 1_8
          ibat(nbat) = irem(j)
        else
          jb = jb + 1_8
          inext(jb) = irem(j)
        end if
      end do

      if (nbat > 0_8) then
        block
          real(r64), allocatable :: txb(:,:), Ab(:,:), Bb(:,:), Cb(:,:)
          integer(8) :: nl, i
          nl = lv%lev(ilev)%n
          allocate(txb(3,nbat), Ab(nbat,nl), Bb(nbat,nl), Cb(nbat,npols))
          do i = 1, nbat
            txb(:,i) = tx(:,ibat(i))
          end do

          if (isimd >= 512_8) then
            call csimd512lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                  lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
          else if (isimd >= 256_8) then
            call csimd256lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                  lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
          else if (isimd >= 128_8) then
            call csimd128lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                  lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
          else
            call lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                  lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
          end if
          call dgemm('N', 'N', nbat, npols, nl, 1.0_r64, Ab, nbat, &
                     REF%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
          do i = 1, nbat
            As(ibat(i),:) = Cb(i,:)
          end do
          call dgemm('N', 'N', nbat, npols, nl, 1.0_r64, Bb, nbat, &
                     REF%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
          do i = 1, nbat
            Ad(ibat(i),:) = Cb(i,:)
          end do

          deallocate(txb, Ab, Bb, Cb)
        end block
      end if

      nrem = jb
      irem(1:nrem) = inext(1:nrem)
    end do

    ms = nrem
    idxs(1:ms) = irem(1:ms)

  end subroutine lap3dsdlpmat_levels_r64

  subroutine adaptive_pass(m, tx, npols, lv, isimd, As, Ad, irem, nrem)
    integer(8),           intent(in)    :: m, npols, isimd
    real(r64),            intent(in)    :: tx(3,m)
    type(patch_levels_t), intent(in)    :: lv
    real(r64),            intent(inout) :: As(m,npols), Ad(m,npols)
    integer(8),           intent(inout) :: irem(m), nrem

    integer(8) :: ntree, n1, j, jj, i, q, il, i0, nl, jb, nb, ic, ig, ngrp, ng
    integer(8) :: cl(64), cs(64), nc, amask(m), keep(m), nkeep, ib(m), inext(m)
    integer(8) :: gmask(NCUTC), key, ngmax
    logical    :: ok

    ntree = lv%nlevel - 1_8
    n1    = REF%lev(1)%n
    ngmax = (4_8**(ntree-1_8))*n1

    key = REF%korder + 1000_8*REF%npols + 1000000_8*REF%order
    if (cutc_key /= key) then
      do i = 1, NCUTC
        CUTC(i)%mask = -1_8
        if (allocated(CUTC(i)%gI)) deallocate(CUTC(i)%gI)
      end do
      ncutn = 0_8;  cutc_key = key
    end if

    if (allocated(WS_Ab)) then
      if (size(WS_Ab,1) < m .or. size(WS_Ab,2) < ngmax .or. &
          size(WS_Cb,2) < npols) deallocate(WS_Ab, WS_Bb, WS_Cb)
    end if
    if (.not. allocated(WS_Ab)) &
      allocate(WS_Ab(m,ngmax), WS_Bb(m,ngmax), WS_Cb(m,npols))

    nkeep = 0_8;  jb = 0_8;  ngrp = 0_8
    do jj = 1, nrem
      j = irem(jj)
      call walk(tx(:,j), ntree, lv, cl, cs, nc, amask(j), ok)
      if (ok) then
        nkeep = nkeep + 1_8;  keep(nkeep) = j
        i = 0_8
        do ig = 1, ngrp
          if (gmask(ig) == amask(j)) then;  i = 1_8;  exit;  end if
        end do
        if (i == 0_8 .and. ngrp < NCUTC) then
          ngrp = ngrp + 1_8;  gmask(ngrp) = amask(j)
        end if
      else
        jb = jb + 1_8;  inext(jb) = j
      end if
    end do

    do ig = 1, ngrp
      nb = 0_8
      do i = 1, nkeep
        if (amask(keep(i)) == gmask(ig)) then;  nb = nb + 1_8;  ib(nb) = keep(i);  end if
      end do
      if (nb == 0_8) cycle

      call walk(tx(:,ib(1)), ntree, lv, cl, cs, nc, key, ok)
      ng = nc*n1
      call cut_gI(gmask(ig), nc, cl, cs, n1, npols, lv, ic)

      block
        real(r64) :: txb(3,nb)
        do i = 1, nb
          txb(:,i) = tx(:,ib(i))
        end do
        do i = 1, nc
          il = cl(i);  q = cs(i);  i0 = (q-1_8)*n1
          if (isimd >= 128_8) then
            call csimd128lap3dsdlpmat_r64(nb, txb, n1, lv%lev(il)%sx(1,i0+1_8), &
                 lv%lev(il)%snx(1,i0+1_8), lv%lev(il)%sw(i0+1_8), &
                 WS_Ab(1:nb,(i-1_8)*n1+1_8:i*n1), WS_Bb(1:nb,(i-1_8)*n1+1_8:i*n1))
          else
            call lap3dsdlpmat_r64(nb, txb, n1, lv%lev(il)%sx(1,i0+1_8), &
                 lv%lev(il)%snx(1,i0+1_8), lv%lev(il)%sw(i0+1_8), &
                 WS_Ab(1:nb,(i-1_8)*n1+1_8:i*n1), WS_Bb(1:nb,(i-1_8)*n1+1_8:i*n1))
          end if
        end do
      end block

      call dgemm('N','N', nb, npols, ng, 1.0_r64, WS_Ab, m, &
                 CUTC(ic)%gI, ng, 0.0_r64, WS_Cb, m)
      do i = 1, nb
        As(ib(i),:) = WS_Cb(i,:)
      end do
      call dgemm('N','N', nb, npols, ng, 1.0_r64, WS_Bb, m, &
                 CUTC(ic)%gI, ng, 0.0_r64, WS_Cb, m)
      do i = 1, nb
        Ad(ib(i),:) = WS_Cb(i,:)
      end do
    end do

    nrem = jb
    irem(1:nrem) = inext(1:nrem)
  end subroutine adaptive_pass

  subroutine cut_gI(mask, nc, cl, cs, n1, npols, lv, ic)
    integer(8),           intent(in)  :: mask, nc, cl(*), cs(*), n1, npols
    type(patch_levels_t), intent(in)  :: lv
    integer(8),           intent(out) :: ic
    integer(8) :: i, k, il, q, i0

    do ic = 1, ncutn
      if (CUTC(ic)%mask == mask) return
    end do
    ncutn = ncutn + 1_8;  ic = ncutn
    CUTC(ic)%mask = mask;  CUTC(ic)%ng = nc*n1
    allocate(CUTC(ic)%gI(nc*n1, npols))
    do i = 1, nc
      il = cl(i);  q = cs(i);  i0 = (q-1_8)*n1
      do k = 1, n1
        CUTC(ic)%gI((i-1_8)*n1+k, :) = REF%lev(il)%interpmat(i0+k, :)
      end do
    end do
  end subroutine cut_gI

  subroutine walk(x, ntree, lv, cl, cs, nc, key, ok)
    real(r64),            intent(in)  :: x(3)
    integer(8),           intent(in)  :: ntree
    type(patch_levels_t), intent(in)  :: lv
    integer(8),           intent(out) :: cl(*), cs(*), nc, key
    logical,              intent(out) :: ok
    integer(8) :: act(64), nxt(64), na, nn, ic, q, b, il, off

    na = 1_8;  act(1) = 1_8;  nc = 0_8;  key = 0_8;  ok = .true.
    do il = 1, ntree
      off = 0_8
      do b = 1, il-1_8
        off = off + 4_8**(b-1_8)
      end do
      nn = 0_8
      do ic = 1, na
        q = act(ic)
        if (sum((x - lv%lev(il)%qpoint(:,q))**2) > lv%lev(il)%qradii2(q)) then
          nc = nc + 1_8;  cl(nc) = il;  cs(nc) = q
          key = ior(key, ishft(1_8, off + q - 1_8))
        else if (il == ntree) then
          ok = .false.;  return
        else
          do b = 1, 4
            nn = nn + 1_8;  nxt(nn) = (q-1_8)*4_8 + b
          end do
        end if
      end do
      if (nn == 0_8) exit
      act(1:nn) = nxt(1:nn);  na = nn
    end do
  end subroutine walk

  subroutine one_level(ilev, m, tx, npols, lv, isimd, As, Ad, irem, nrem)
    integer(8),           intent(in)    :: ilev, m, npols, isimd
    real(r64),            intent(in)    :: tx(3,m)
    type(patch_levels_t), intent(in)    :: lv
    real(r64),            intent(inout) :: As(m,npols), Ad(m,npols)
    integer(8),           intent(inout) :: irem(m), nrem

    integer(8) :: j, k, jb, nbat, ibat(m), inext(m), nl, i
    logical    :: outside

    nbat = 0_8;  jb = 0_8
    do j = 1, nrem
      outside = .true.
      do k = 1, lv%lev(ilev)%nsub
        if (sum((lv%lev(ilev)%qpoint(:,k) - tx(:,irem(j)))**2) &
            < lv%lev(ilev)%qradii2(k)) then
          outside = .false.;  exit
        end if
      end do
      if (outside) then
        nbat = nbat + 1_8;  ibat(nbat) = irem(j)
      else
        jb = jb + 1_8;  inext(jb) = irem(j)
      end if
    end do

    if (nbat > 0_8) then
      block
        real(r64), allocatable :: txb(:,:), Ab(:,:), Bb(:,:), Cb(:,:)
        nl = lv%lev(ilev)%n
        allocate(txb(3,nbat), Ab(nbat,nl), Bb(nbat,nl), Cb(nbat,npols))
        do i = 1, nbat
          txb(:,i) = tx(:,ibat(i))
        end do
        if (isimd >= 512_8) then
          call csimd512lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
        else if (isimd >= 256_8) then
          call csimd256lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
        else if (isimd >= 128_8) then
          call csimd128lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
        else
          call lap3dsdlpmat_r64(nbat, txb, nl, lv%lev(ilev)%sx, &
                                lv%lev(ilev)%snx, lv%lev(ilev)%sw, Ab, Bb)
        end if
        call dgemm('N','N', nbat, npols, nl, 1.0_r64, Ab, nbat, &
                   REF%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
        do i = 1, nbat
          As(ibat(i),:) = Cb(i,:)
        end do
        call dgemm('N','N', nbat, npols, nl, 1.0_r64, Bb, nbat, &
                   REF%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
        do i = 1, nbat
          Ad(ibat(i),:) = Cb(i,:)
        end do
        deallocate(txb, Ab, Bb, Cb)
      end block
    end if

    nrem = jb
    irem(1:nrem) = inext(1:nrem)
  end subroutine one_level

end module patch_refine_mod
