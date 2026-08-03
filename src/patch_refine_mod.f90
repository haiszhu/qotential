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
            lap3dsdlpmat_levels_r64

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
    real(r64), allocatable :: uv(:,:)        ! (2,n)  reference coordinates
    real(r64), allocatable :: wt(:)          ! (n)    reference weights
    real(r64), allocatable :: sx(:,:)        ! (3,n)  surface points
    real(r64), allocatable :: snx(:,:)       ! (3,n)  unit normals
    real(r64), allocatable :: sw(:)          ! (n)    area weights
    real(r64), allocatable :: interpmat(:,:) ! (n,npols) from the base patch
    real(r64), allocatable :: qpoint(:,:)    ! (3,nsub) weighted centroid
    real(r64), allocatable :: qradii2(:)     ! (nsub)   squared near radius
  end type patch_level_t

  type :: patch_levels_t
    integer(8) :: nlevel = 0_8
    type(patch_level_t), allocatable :: lev(:)   ! 1=ff, 2=f, 3=c, 4=cc
  end type patch_levels_t

contains

  subroutine build_patch_levels_r64(korder, npols, sx, rts, rps, order, &
                                    dist, nlevel, lv)
    integer(8),           intent(in)  :: korder, npols, order, nlevel
    real(r64),            intent(in)  :: sx(3,npols), rts(3,npols), rps(3,npols)
    real(r64),            intent(in)  :: dist
    type(patch_levels_t), intent(out) :: lv

    real(r64)  :: umatr(npols,npols), vmatr(npols,npols)
    real(r64)  :: srcT(npols,9), dist2
    integer(8) :: ilev, n1, nprev, n, nsub, ncc

    dist2 = dist*dist

    srcT(:,1:3) = transpose(sx)
    srcT(:,4:6) = transpose(rts)
    srcT(:,7:9) = transpose(rps)

    umatr = 0.0_r64;  vmatr = 0.0_r64
    call koorn_vals2coefs_coefs2vals(korder, npols, umatr, vmatr)

    n1  = order*(order+1_8)/2_8        ! nodes per sub-patch, VR order-1
    ncc = size(CC_SUBIDX, kind=8) + 4_8
    lv%nlevel = nlevel
    allocate(lv%lev(nlevel))

    do ilev = 1, nlevel

      if (ilev == 1_8) then
        n = n1;  nsub = 1_8
        allocate(lv%lev(1)%uv(2,n), lv%lev(1)%wt(n))
        call get_vioreanu_nodes(order-1_8, n, lv%lev(1)%uv)
        call get_vioreanu_wts  (order-1_8, n, lv%lev(1)%wt)

      else if (ilev == nlevel .and. nlevel == 4_8) then
        n = ncc*n1;  nsub = ncc
        allocate(lv%lev(ilev)%uv(2,n), lv%lev(ilev)%wt(n))
        block
          real(r64), allocatable :: uvt(:,:), wtt(:)
          real(r64)  :: uv0(2)
          integer(8) :: np, q, b, i, isrc, idst
          np = lv%lev(ilev-1)%n
          allocate(uvt(2,4*np), wtt(4*np))
          do i = 1, np
            uv0 = 2.0_r64*lv%lev(ilev-1)%uv(:,i) - 1.0_r64
            uvt(:,      i) =  uv0
            uvt(:,   np+i) = -uv0
            uvt(1, 2*np+i) =  uv0(1)
            uvt(2, 2*np+i) =  uv0(2) + 2.0_r64
            uvt(1, 3*np+i) =  uv0(1) + 2.0_r64
            uvt(2, 3*np+i) =  uv0(2)
            do b = 0, 3
              wtt(b*np+i) = 0.25_r64*lv%lev(ilev-1)%wt(i)
            end do
          end do
          uvt = ((uvt - 1.0_r64)*0.5_r64 + 1.0_r64)*0.5_r64
          do q = 1, size(CC_SUBIDX, kind=8)
            do i = 1, n1
              isrc = (CC_SUBIDX(q)-1_8)*n1 + i
              idst = (q-1_8)*n1 + i
              lv%lev(ilev)%uv(:,idst) = uvt(:,isrc)
              lv%lev(ilev)%wt(idst)   = wtt(isrc)
            end do
          end do
          do q = 1, 4
            do i = 1, n1
              isrc = (q-1_8)*n1 + i
              idst = (size(CC_SUBIDX, kind=8) + q - 1_8)*n1 + i
              lv%lev(ilev)%uv(:,idst) = &
                   (lv%lev(2)%uv(:,isrc) - 0.5_r64)*CC_EXTRA_SCALE + CC_EXTRA_OFFSET
              lv%lev(ilev)%wt(idst)   = CC_EXTRA_SCALE**2 * lv%lev(2)%wt(isrc)
            end do
          end do
          deallocate(uvt, wtt)
        end block

      else
        nprev = lv%lev(ilev-1)%n
        n = 4_8*nprev;  nsub = 4_8*lv%lev(ilev-1)%nsub
        allocate(lv%lev(ilev)%uv(2,n), lv%lev(ilev)%wt(n))
        block
          real(r64)  :: uv0(2)
          integer(8) :: i, b
          do i = 1, nprev
            uv0 = 2.0_r64*lv%lev(ilev-1)%uv(:,i) - 1.0_r64
            lv%lev(ilev)%uv(:,         i) =  uv0
            lv%lev(ilev)%uv(:,  nprev+i) = -uv0
            lv%lev(ilev)%uv(1,2*nprev+i) =  uv0(1)
            lv%lev(ilev)%uv(2,2*nprev+i) =  uv0(2) + 2.0_r64
            lv%lev(ilev)%uv(1,3*nprev+i) =  uv0(1) + 2.0_r64
            lv%lev(ilev)%uv(2,3*nprev+i) =  uv0(2)
            do b = 0, 3
              lv%lev(ilev)%wt(b*nprev+i) = 0.25_r64*lv%lev(ilev-1)%wt(i)
            end do
          end do
          lv%lev(ilev)%uv = ((lv%lev(ilev)%uv - 1.0_r64)*0.5_r64 + 1.0_r64)*0.5_r64
        end block
      end if

      lv%lev(ilev)%n    = n
      lv%lev(ilev)%nsub = nsub

      allocate(lv%lev(ilev)%sx(3,n), lv%lev(ilev)%snx(3,n), lv%lev(ilev)%sw(n))
      allocate(lv%lev(ilev)%interpmat(n,npols))
      allocate(lv%lev(ilev)%qpoint(3,nsub), lv%lev(ilev)%qradii2(nsub))

      block
        real(r64), allocatable :: pols(:,:), geo(:,:)
        real(r64)  :: sp, swsum, qp(3), d2, rmax2
        integer(8) :: i, k, i0

        allocate(pols(n,npols), geo(n,9))
        call koorn_pols_batch_r64(n, lv%lev(ilev)%uv, korder, npols, pols)

        call dgemm('N', 'N', n, npols, npols, 1.0_r64, pols, n, &
                   umatr, npols, 0.0_r64, lv%lev(ilev)%interpmat, n)
        call dgemm('N', 'N', n, 9_8, npols, 1.0_r64, lv%lev(ilev)%interpmat, n, &
                   srcT, npols, 0.0_r64, geo, n)

        lv%lev(ilev)%sx = transpose(geo(:,1:3))
        do i = 1, n
          lv%lev(ilev)%snx(1,i) = geo(i,8)*geo(i,6) - geo(i,9)*geo(i,5)
          lv%lev(ilev)%snx(2,i) = geo(i,9)*geo(i,4) - geo(i,7)*geo(i,6)
          lv%lev(ilev)%snx(3,i) = geo(i,7)*geo(i,5) - geo(i,8)*geo(i,4)
          sp = sqrt(dot_product(lv%lev(ilev)%snx(:,i), lv%lev(ilev)%snx(:,i)))
          lv%lev(ilev)%snx(:,i) = lv%lev(ilev)%snx(:,i)/sp
          lv%lev(ilev)%sw(i)    = sp*lv%lev(ilev)%wt(i)
        end do

        do k = 1, nsub
          i0    = (k-1_8)*n1
          swsum = sum(lv%lev(ilev)%sw(i0+1:i0+n1))
          qp(1) = dot_product(lv%lev(ilev)%sx(1,i0+1:i0+n1), &
                              lv%lev(ilev)%sw(i0+1:i0+n1))/swsum
          qp(2) = dot_product(lv%lev(ilev)%sx(2,i0+1:i0+n1), &
                              lv%lev(ilev)%sw(i0+1:i0+n1))/swsum
          qp(3) = dot_product(lv%lev(ilev)%sx(3,i0+1:i0+n1), &
                              lv%lev(ilev)%sw(i0+1:i0+n1))/swsum
          rmax2 = 0.0_r64
          do i = i0+1, i0+n1
            d2 = sum((lv%lev(ilev)%sx(:,i) - qp)**2)
            rmax2 = max(rmax2, d2)
          end do
          lv%lev(ilev)%qpoint(:,k) = qp
          lv%lev(ilev)%qradii2(k)  = dist2*max(rmax2, swsum)
        end do

        deallocate(pols, geo)
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
                     lv%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
          do i = 1, nbat
            As(ibat(i),:) = Cb(i,:)
          end do
          call dgemm('N', 'N', nbat, npols, nl, 1.0_r64, Bb, nbat, &
                     lv%lev(ilev)%interpmat, nl, 0.0_r64, Cb, nbat)
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

end module patch_refine_mod
