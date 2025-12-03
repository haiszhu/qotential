      ! -----------------------------------------------------------------
      subroutine Lap3dDLP_closepaneladp_vr_guru(m,tx,
     1                                          n,sx,sw,snx,iside,A)
      implicit none
      integer *8, intent(in) :: m, n, iside
      real *8, intent(in) :: tx(3,m), sx(3,n), sw(n), snx(3,n)
      real *8, intent(inout) :: A(m,n)
      ! 
      integer *8 :: sbdnp, len, nbd, nquad, isimd
      integer *8 :: hdim, nterms
      real *8 :: umatr(n,n), alpha, sxc(3), r_vert(3,3)
      real *8 :: Mmatrix(4*n,4*n), Omega(4*n,m)
      ! 
      real *8, allocatable :: sxbd(:,:), stangbd(:,:)
      real *8, allocatable :: swbd(:), sspbd(:)
      real *8, allocatable :: tgl(:), wgl(:), w_bclag(:)
      real *8, allocatable :: Dgl(:,:), Agl(:,:)
      real *8, allocatable :: Legmat(:,:), bclagmatlr(:,:)
      real *8, allocatable :: onm0(:,:,:), onm1(:,:,:)
      real *8, allocatable :: onm2(:,:,:), onm3(:,:,:)
      complex *16, allocatable :: Fbd(:,:)
      complex *16, allocatable :: Fxbd(:,:), Fybd(:,:), Fzbd(:,:)
      
      ! 
      nterms = int(sqrt(dble(1+8*n)) + 0.5d0, kind=8)
      nterms = (nterms - 1)/2
      alpha = 0
      isimd = 0
      nquad = nterms+8
      len = 3
      sbdnp = 3*len
      nbd = sbdnp*nquad
      ! 
      allocate(sxbd(3,nbd), stangbd(3,nbd))
      allocate(swbd(nbd), sspbd(nbd))
      allocate(tgl(nquad), wgl(nquad), w_bclag(nquad))
      allocate(Dgl(nquad,nquad), Agl(nquad,nquad))
      allocate(Legmat(nquad,nquad), bclagmatlr(nquad,2))
      allocate(onm0(4,n,nbd), onm1(4,n,nbd))
      allocate(onm2(4,n,nbd), onm3(4,n,nbd))
      allocate(Fbd(nbd,(nterms+1)*(nterms+1)))
      allocate(Fxbd(nbd,(nterms+1)*(nterms+1)))
      allocate(Fybd(nbd,(nterms+1)*(nterms+1)))
      allocate(Fzbd(nbd,(nterms+1)*(nterms+1)))
      call Lap3dDLP_closepaneladp_vr(m,tx,nterms,
     1                      n,sx,sw,snx,umatr,iside,hdim,
     2                      alpha,sxc,r_vert,
     3                      sbdnp,len,nbd,nquad,sxbd,stangbd,swbd,sspbd,
     4                      tgl,wgl,w_bclag,Dgl,Agl,Legmat,bclagmatlr,
     5                      onm0,onm1,onm2,onm3,Fbd,Fxbd,Fybd,Fzbd,
     6                      Mmatrix,isimd,A,Omega)

      deallocate(sxbd, stangbd)
      deallocate(swbd, sspbd)
      deallocate(tgl, wgl, w_bclag)
      deallocate(Dgl, Agl)
      deallocate(Legmat, bclagmatlr)
      deallocate(onm0, onm1, onm2, onm3)
      deallocate(Fbd, Fxbd, Fybd, Fzbd)
      end subroutine Lap3dDLP_closepaneladp_vr_guru
      
      ! simplest version with no optimization on parameters
      ! the interface is also for personal use...
      subroutine Lap3dDLP_closepaneladp_vr(m,tx,nterms,
     1                      n,sx,sw,snx,umatr,iside,hdim,
     2                      alpha,sxc,r_vert,
     3                      sbdnp,len,nbd,nquad,sxbd,stangbd,swbd,sspbd,
     4                      tgl,wgl,w_bclag,Dgl,Agl,Legmat,bclagmatlr,
     5                      onm0,onm1,onm2,onm3,Fbd,Fxbd,Fybd,Fzbd,
     6                      Mmatrix,isimd,A,Omega)
      implicit none
      ! Declarations from original code (assumed correct)
      integer *8, intent(in) :: m, n, nterms, iside
      integer *8, intent(in) :: sbdnp, len, nbd, nquad, isimd
      real *8, intent(in) :: tx(3,m), sx(3,n), sw(n), snx(3,n)
      integer *8, intent(inout) :: hdim
      real *8, intent(inout) :: umatr(n,n), alpha, sxc(3), r_vert(3,3)
      real *8, intent(inout) :: sxbd(3,nbd), stangbd(3,nbd)
      real *8, intent(inout) :: swbd(nbd), sspbd(nbd)
      real *8, intent(inout) :: tgl(nquad), wgl(nquad), w_bclag(nquad)
      real *8, intent(inout) :: Dgl(nquad,nquad), Agl(nquad,nquad)
      real *8, intent(inout) :: Legmat(nquad,nquad), bclagmatlr(nquad,2)
      real *8, intent(inout) :: onm0(4,n,nbd), onm1(4,n,nbd)
      real *8, intent(inout) :: onm2(4,n,nbd), onm3(4,n,nbd)
      complex *16, intent(inout) :: Fbd(nbd,(nterms+1)*(nterms+1))
      complex *16, intent(inout) :: Fxbd(nbd,(nterms+1)*(nterms+1))
      complex *16, intent(inout) :: Fybd(nbd,(nterms+1)*(nterms+1))
      complex *16, intent(inout) :: Fzbd(nbd,(nterms+1)*(nterms+1))
      real *8, intent(inout) :: Mmatrix(4*n,4*n), A(m,n), Omega(4*n,m)

      real *8 :: o_offset, pi
      integer *8 :: itype, k, j, ier, korder
      real *8 :: v(nquad,nquad), vp(nquad,nquad), vpp(nquad,nquad)
      real *8 :: sx_mean(3)
      real *8 :: txnew0(3,m), sxnew0(3,n), snxnew0(3,n)
      real *8 :: vmatr(n,n)
      real *8 :: lsstep, tpan(sbdnp+1)
      real *8 :: r1(3), r2(3), r3(3), nc(3), temp
      real *8 :: txnew(3,m), sxnew(3,n), snxnew(3,n)
      real *8 :: sx3min, sx3max
      integer *8 :: idxvec(n), tmpidx, tmpidx2, ij
      real *8 :: F(nbd,n), F1(nbd,n), F2(nbd,n), F3(nbd,n)
      integer *8 :: morder
      real *8, allocatable :: m_all_adp(:,:,:)
      real *8 :: tglr,tgll,denoml,denomr,xdiffl,xdiffr,templ,tempr
      real *8 :: rho, rat1(2,nquad)
      integer *8 :: ell, idx_ell_start, idx_ell_end
      real *8 :: DglT(nquad,nquad), r_ell(3,nquad), rp(3,nquad),w(nquad)
      real *8 :: rlr(3,2), rl(3), rr(3), pan_len
      real *8 :: r0j(3), x0, y0, z0
      integer *8 :: rfc
      real *8 :: troot_real, xroot_real, yroot_real, zroot_real
      real *8 :: r_root(3), sqn_dist
      integer *8 :: lenj, lenl, lenr, nquad_up
      real *8, allocatable :: r_up(:,:), rp_up(:,:), Br(:,:), w_up(:)
      integer *8 :: flag
      real *8, allocatable :: mk_j_up(:,:), nk_j_up(:,:)
      real *8 :: mk_j(nquad,nterms+2), nk_j(nquad,nterms+2)
      integer *8 :: ijIdx(2,n), idx, idx_k_end, dim1
      real *8 :: omega0(n,m),omega1(n,m),omega2(n,m),omega3(n,m)
      complex *16 :: Fc(n,(nterms+1)*(nterms+1))
      complex *16 :: Fx_c(n,(nterms+1)*(nterms+1))
      complex *16 :: Fy_c(n,(nterms+1)*(nterms+1))
      complex *16 :: Fz_c(n,(nterms+1)*(nterms+1))
      integer *8 :: info
      real *8 :: F0sx(n,n), F1sx(n,n), F2sx(n,n), F3sx(n,n)
      real *8 :: MmatrixT(4*n,4*n), Atmp(4*n,m)

      pi = 4.0d0*atan(1.0d0)
      o_offset = 0.063d0
      hdim = nterms*(nterms+1)/2

      ! 
      tgl = 0.0d0
      wgl = 0.0d0
      Dgl = 0.0d0
      w_bclag = 0.0d0
      Legmat = 0.0d0
      Agl = 0.0d0
      v = 0.0d0
      vp = 0.0d0
      vpp = 0.0d0
      itype = int(2, KIND=8)
      call legeexps2(itype,nquad,tgl,Legmat,v,wgl,vp,vpp);
      call bclaginterpweights(nquad,tgl,w_bclag);
      Dgl = matmul(vp,Legmat)
      Agl = 0.0d0
      Agl(1,:) = 1.0d0
      do k=2,nquad
         do j=1,nquad
            Agl(k,j) = Agl(k-1,j)*tgl(j)
         end do
      end do
      ier = 0

      ! 
      sx_mean = 0.0d0
      do k=1,n
          sx_mean = sx_mean + sx(:,k)
      end do
      sx_mean = sx_mean / dble(n)
      alpha = 1.25d0 / sqrt(sum(sw)*2.0d0)
      do k=1,m
         txnew0(1,k) = alpha*(tx(1,k)-sx_mean(1))
         txnew0(2,k) = alpha*(tx(2,k)-sx_mean(2))
         txnew0(3,k) = alpha*(tx(3,k)-sx_mean(3))
      end do
      do k=1,n
         sxnew0(1,k) = alpha*(sx(1,k)-sx_mean(1))
         sxnew0(2,k) = alpha*(sx(2,k)-sx_mean(2))
         sxnew0(3,k) = alpha*(sx(3,k)-sx_mean(3))
      end do
      snxnew0 = snx

      ! 
      korder = nterms-1
      vmatr = 0.0d0
      umatr = 0.0d0
      call koorn_vals2coefs_coefs2vals(korder,n,umatr,vmatr)

      ! 
      lsstep = 2*pi/sbdnp
      tpan = 0.0d0
      do k=1,sbdnp+1
        tpan(k) = (k-1)*lsstep
      enddo
      sxbd = 0.0d0
      swbd = 0.0d0
      stangbd = 0.0d0
      sspbd = 0.0d0
      r_vert = 0.0d0
      call line3quadr_3dline(sxnew0, korder, n, umatr, 
     1                  nquad, tgl, wgl, Dgl, sbdnp, tpan, nbd, 
     2                  sxbd, swbd, stangbd, sspbd, r_vert)

      ! 
      r1 = r_vert(:,1)
      r2 = r_vert(:,2)
      r3 = r_vert(:,3)
      nc(1) = (r2(2)-r1(2))*(r3(3)-r1(3)) - (r2(3)-r1(3))*(r3(2)-r1(2))
      nc(2) = (r2(3)-r1(3))*(r3(1)-r1(1)) - (r2(1)-r1(1))*(r3(3)-r1(3))
      nc(3) = (r2(1)-r1(1))*(r3(2)-r1(2)) - (r2(2)-r1(2))*(r3(1)-r1(1))
      temp = sqrt(nc(1)**2 + nc(2)**2 + nc(3)**2)
      if (temp .gt. 0.0d0) nc = nc/temp
      sxnew = 0.0d0
      call transformedsx(r_vert,n,sxnew0,sxnew);
      if (iside .eq. 0) then 
         sx3min = max(abs(minval(sxnew(3,:))),o_offset)
         do j=1,3
            do k=1,3
               r_vert(j,k) = r_vert(j,k) - 4.25d0*sx3min*nc(j)
            end do
         end do
      else if (iside .eq. 1) then 
         sx3max = max(abs(maxval(sxnew(3,:))),o_offset)
         do j=1,3
            do k=1,3
               r_vert(j,k) = r_vert(j,k) + 4.25d0*sx3max*nc(j)
            end do
         end do
      end if

      snxnew0 = snx
      call transformedsxsnxtx(r_vert,n,sxnew0,snxnew0,m,txnew0)
      sxnew = sxnew0
      snxnew = snxnew0
      txnew = txnew0

      call line3quadr_3dline(sxnew, korder, n, umatr, nquad, tgl, wgl, 
     1      Dgl, sbdnp, tpan, nbd, sxbd, swbd, stangbd, sspbd, r_vert)

      idxvec = 0
      tmpidx = 0
      tmpidx2 = 0
      do ij = 0,nterms
        do k = -ij,ij
          tmpidx2 = tmpidx2 + 1
          if ((ij.gt.0).and.(k.gt.0)) then
            tmpidx = tmpidx + 1
            idxvec(tmpidx) = tmpidx2
          endif
        enddo
      enddo

      ! 
      Omega = 0.0d0
      Fbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fxbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fybd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fzbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      call l3dtavecevalmatf(sxbd,nbd,nterms,Fbd,Fxbd,Fybd,Fzbd,ier)
      do k = 1,n
        F(:,k) = dble(Fbd(:,idxvec(k)))
        F1(:,k) = dble(Fxbd(:,idxvec(k)))
        F2(:,k) = dble(Fybd(:,idxvec(k)))
        F3(:,k) = dble(Fzbd(:,idxvec(k)))
      enddo
      
      ! 
      onm0 = 0.0d0
      onm1 = 0.0d0
      onm2 = 0.0d0
      onm3 = 0.0d0
      call omeganm_precompff(nbd,n,F,F1,F2,F3,sxbd,stangbd,swbd,
     1                      onm0,onm1,onm2,onm3)
      
      ! 
      morder = nterms+2
      allocate( m_all_adp(morder,nbd,m))
      m_all_adp = 0.0d0

      tglr =  1.0d0
      tgll = -1.0d0
      bclagmatlr = 0.0d0
      denoml = 0.0d0
      denomr = 0.0d0
      do j=1,nquad
        xdiffl = tgll-tgl(j)
        xdiffr = tglr-tgl(j)
        templ = w_bclag(j)/xdiffl
        tempr = w_bclag(j)/xdiffr
        bclagmatlr(j,1) = templ
        bclagmatlr(j,2) = tempr
        denoml = denoml + templ
        denomr = denomr + tempr
      enddo
      denoml = 1.0D0/denoml
      denomr = 1.0D0/denomr
      bclagmatlr(:,1) = bclagmatlr(:,1) * denoml
      bclagmatlr(:,2) = bclagmatlr(:,2) * denomr

      ! 
      rho = 4.0d0**(16.0d0/dble(nquad))
      rat1 = 0.0d0
      call legendrescalarf_init(nquad-1,rat1)
      DglT = transpose(Dgl)
      do ell = 1,sbdnp
        idx_ell_start = (ell-1)*nquad + 1
        idx_ell_end = ell*nquad
        r_ell = sxbd(:, idx_ell_start:idx_ell_end)
        rp = matmul(r_ell,DglT)
        w = sqrt(rp(1,:)*rp(1,:)+rp(2,:)*rp(2,:)+rp(3,:)*rp(3,:))*wgl
        rlr = matmul(r_ell,bclagmatlr)
        rl = rlr(:,1)
        rr = rlr(:,2)
        pan_len = sum(w)
        ! 
        do j = 1,m
          r0j = txnew(:, j)
          x0 = r0j(1)
          y0 = r0j(2)
          z0 = r0j(3)
          ! 
          rfc = int(0, KIND=8) ! initialize rfc to false
          troot_real = 0.0d0
          xroot_real = 0.0d0
          yroot_real = 0.0d0
          zroot_real = 0.0d0
          call line3nearrootf_real(tgl, Agl, Legmat, wgl, 
     1   r_ell(1,:), r_ell(2,:), r_ell(3,:), nquad, x0, y0, z0, 
     2   rho, rfc, troot_real, xroot_real, yroot_real, zroot_real, rat1)

          if (rfc>0) then
            r_root(1) = xroot_real
            r_root(2) = yroot_real
            r_root(3) = zroot_real
            sqn_dist = 0.0d0
            lenj = int(0, KIND=8)
            lenl = int(0, KIND=8)
            lenr = int(0, KIND=8)
            ! 
            call lenofpantendforsaf(nquad,tgl,wgl,r_ell,rr,rl,r0j,
     1              troot_real,r_root,w_bclag,Dgl,sqn_dist,
     2              pan_len,lenj,lenl,lenr,rp,w)
            ! 
            nquad_up = nquad*(lenj-1)
            allocate(r_up(3,nquad_up),rp_up(3,nquad_up))
            allocate(Br(nquad,nquad_up),w_up(nquad_up))
            r_up = 0.0d0
            rp_up = 0.0d0
            Br = 0.0d0
            w_up = 0.0d0
            call line3adaptivernearr0forsaf(troot_real, nquad, tgl, wgl, 
     1              r_ell, w_bclag, lenj, lenl, lenr, rp, w, nquad_up, 
     2              r_up, rp_up, Br, w_up)
            ! 
            flag = int(0, KIND=8)
            allocate(mk_j_up(nquad_up,morder))
            allocate(nk_j_up(nquad_up,morder))
            mk_j_up = 0.0d0
            nk_j_up = 0.0d0
            call momentsad_vr(r0j, nquad_up, r_up, nterms+1, flag, 
     1              nk_j_up, mk_j_up)
            mk_j = matmul(Br, mk_j_up)

            deallocate(r_up,rp_up)
            deallocate(Br,w_up)
            deallocate(mk_j_up)
            deallocate(nk_j_up)
          else
            ! 
            flag = int(0, KIND=8)
            mk_j = 0.0d0 
            nk_j = 0.0d0
            call momentsad_vr(r0j, nquad, r_ell, nterms+1, flag, 
     1              nk_j, mk_j)
          end if  
          m_all_adp(1:morder,idx_ell_start:idx_ell_end,j) = 
     1              transpose(mk_j)
        enddo
      enddo

      ! 
      ijIdx = int(0, KIND=8)
      idx = int(1, KIND=8)
      do k=1,nterms
        idx_k_end = idx + k - 1
        ijIdx(1, idx:idx_k_end) = ijIdx(1, idx:idx_k_end) + (k - 1)
        idx = idx + k
      enddo
      omega0 = 0.0d0
      omega1 = 0.0d0
      omega2 = 0.0d0
      omega3 = 0.0d0
      call omegaall0123vrnewf(m, txnew, 
     1         m_all_adp, nbd, nterms, morder, onm0, onm1, onm2, onm3, 
     2         hdim, ijIdx, omega0, omega1, omega2, omega3)
      do k=1,m
         Omega(1:n,k) = omega0(1:n,k)
         Omega(n+1:2*n,k) = -omega1(1:n,k)
         Omega(2*n+1:3*n,k) = -omega2(1:n,k)
         Omega(3*n+1:4*n,k) = -omega3(1:n,k)
      end do

      Fc = (0.0d0,0.0d0)
      Fx_c = (0.0d0,0.0d0)
      Fy_c = (0.0d0,0.0d0)
      Fz_c = (0.0d0,0.0d0)
      info = 0
      call l3dtavecevalmatf(sxnew, n, nterms, Fc, Fx_c, Fy_c, Fz_c,info)

      ! 
      F0sx = 0.0d0
      do k=1,hdim
        F1sx(:,k) = dble(Fx_c(:,idxvec(k)))
        F2sx(:,k) = dble(Fy_c(:,idxvec(k)))
        F3sx(:,k) = dble(Fz_c(:,idxvec(k)))
      end do
      Mmatrix = 0.0d0
      Mmatrix(1:n,(hdim+1):(2*hdim)) = -F1sx
      Mmatrix(1:n,(2*hdim+1):(3*hdim)) = -F2sx
      Mmatrix(1:n,(3*hdim+1):(4*hdim)) = -F3sx
      Mmatrix((n+1):(2*n),1:hdim) = F1sx
      Mmatrix((n+1):(2*n),(2*hdim+1):(3*hdim)) = -F3sx
      Mmatrix((n+1):(2*n),(3*hdim+1):(4*hdim)) = F2sx
      Mmatrix((2*n+1):(3*n),1:hdim) = F2sx
      Mmatrix((2*n+1):(3*n),(hdim+1):(2*hdim)) = F3sx
      Mmatrix((2*n+1):(3*n),(3*hdim+1):(4*hdim)) = -F1sx
      Mmatrix((3*n+1):(4*n),1:hdim) = F3sx
      Mmatrix((3*n+1):(4*n),(hdim+1):(2*hdim)) = -F2sx
      Mmatrix((3*n+1):(4*n),(2*hdim+1):(3*hdim)) = F1sx
      
      MmatrixT = transpose(Mmatrix)
      Atmp = Omega
      call dgausselimvec(4*n,MmatrixT,m,Atmp,info)
      A = transpose(Atmp(1:n,:))

      deallocate( m_all_adp )

      end subroutine Lap3dDLP_closepaneladp_vr
      
      subroutine omeganm_precompff( nbd, hdim, F, F1, F2, F3,
     1              sxbd, stangbd, swbd, onm0, onm1, onm2, onm3)
      implicit none
      integer *8, intent(in) :: nbd, hdim
      ! 
      real *8, intent(in) :: F(nbd,hdim)
      real *8, intent(in) :: F1(nbd,hdim), F2(nbd,hdim), F3(nbd,hdim)
      real *8, intent(in) :: sxbd(3,nbd), stangbd(3,nbd), swbd(nbd)
      real *8, intent(inout) :: onm0(4,hdim,nbd), onm1(4,hdim,nbd)
      real *8, intent(inout) :: onm2(4,hdim,nbd), onm3(4,hdim,nbd)

      integer *8 :: i, k, j
      real *8 :: dx(nbd), dy(nbd), dz(nbd)
      real *8 :: f1t, f2t, f3t
      real *8 :: sx1, sx2, sx3
      real *8 :: tx, ty, tz
      real *8 :: sx1f1, sx2f2, sx3f3, mr_dot_f
      real *8 :: q01(4), q02(4), q03(4)
      real *8 :: q11(4), q12(4), q13(4)
      real *8 :: q21(4), q22(4), q23(4)
      real *8 :: q31(4), q32(4), q33(4)
      real *8 :: c1, c2, c3, c4, c5, c6
      real *8 :: acoef, bcoef, ccoef

      ! --------
      do k = 1, nbd
         dx(k) = stangbd(1,k)*swbd(k)
         dy(k) = stangbd(2,k)*swbd(k)
         dz(k) = stangbd(3,k)*swbd(k)
      end do

      ! --------
      onm0 = 0.0d0
      onm1 = 0.0d0
      onm2 = 0.0d0
      onm3 = 0.0d0

      ! --------
      do k = 1, nbd
         sx1 = sxbd(1,k)
         sx2 = sxbd(2,k)
         sx3 = sxbd(3,k)
         tx  = dx(k)
         ty  = dy(k)
         tz  = dz(k)

         !
         c1 = tx*sx3
         c2 = tx*sx2
         c3 = ty*sx1
         c4 = ty*sx3
         c5 = tz*sx2
         c6 = tz*sx1

         ! 
         acoef = -c4 + c5
         bcoef =  c1 - c6
         ccoef = -c2 + c3

         do i = 1, hdim
            f1t = F1(k,i)
            f2t = F2(k,i)
            f3t = F3(k,i)

            ! 
            sx1f1 = sx1*f1t
            sx2f2 = sx2*f2t
            sx3f3 = sx3*f3t
            mr_dot_f = -(sx1f1 + sx2f2 + sx3f3)

            ! ===== qnm_0_* =====
            q01(1) = -sx2*f3t + sx3*f2t
            q01(2) =  0.0d0
            q01(3) =  f3t
            q01(4) = -f2t

            q02(1) = -sx3*f1t + sx1*f3t
            q02(2) = -f3t
            q02(3) =  0.0d0
            q02(4) =  f1t

            q03(1) = -sx1*f2t + sx2*f1t
            q03(2) =  f2t
            q03(3) = -f1t
            q03(4) =  0.0d0

            ! ===== qnm_1_* =====
            q11(1) =  mr_dot_f + 2.0d0*sx1*f1t
            q11(2) = -f1t
            q11(3) =  f2t
            q11(4) =  f3t

            q12(1) =  sx1*f2t + sx2*f1t
            q12(2) = -f2t
            q12(3) = -f1t
            q12(4) =  0.0d0

            q13(1) =  sx1*f3t + sx3*f1t
            q13(2) = -f3t
            q13(3) =  0.0d0
            q13(4) = -f1t

            ! ===== qnm_2_* =====
            q21(1) =  sx2*f1t + sx1*f2t
            q21(2) = -f2t
            q21(3) = -f1t
            q21(4) =  0.0d0

            q22(1) =  mr_dot_f + 2.0d0*sx2*f2t
            q22(2) =  f1t
            q22(3) = -f2t
            q22(4) =  f3t

            q23(1) =  sx2*f3t + sx3*f2t
            q23(2) =  0.0d0
            q23(3) = -f3t
            q23(4) = -f2t

            ! ===== qnm_3_* =====
            q31(1) =  sx3*f1t + sx1*f3t
            q31(2) = -f3t
            q31(3) =  0.0d0
            q31(4) = -f1t

            q32(1) =  sx3*f2t + sx2*f3t
            q32(2) =  0.0d0
            q32(3) = -f3t
            q32(4) = -f2t

            q33(1) =  mr_dot_f + 2.0d0*sx3*f3t
            q33(2) =  f1t
            q33(3) =  f2t
            q33(4) = -f3t

            ! ==========
            do j = 1, 4
               onm0(j,i,k) = acoef*q01(j) + bcoef*q02(j) + ccoef*q03(j)
               onm1(j,i,k) = acoef*q11(j) + bcoef*q12(j) + ccoef*q13(j)
               onm2(j,i,k) = acoef*q21(j) + bcoef*q22(j) + ccoef*q23(j)
               onm3(j,i,k) = acoef*q31(j) + bcoef*q32(j) + ccoef*q33(j)
            end do

         end do
      end do

      end subroutine omeganm_precompff


      subroutine lenofpantendforsaf(p, tgl, wgl, r, rr, rl, r0, t_root,
     1   r_root, w_bclag, D, sqn_dist2, pan_len, len, lenl, lenr, rp, w)

      implicit none
      integer *8, intent(in) :: p
      real *8, intent(in) :: tgl(p), wgl(p), w_bclag(p), r(3,p)
      real *8, intent(in) :: rr(3), rl(3)
      real *8, intent(in) :: r0(3), r_root(3), t_root, D(p,p)
      real *8, intent(inout) :: sqn_dist2, pan_len, rp(3,p), w(p)
      integer *8, intent(inout) :: len, lenl, lenr

      real *8, parameter :: factor        = 3.0D0
      real *8, parameter :: factor_sq     = factor*factor
      real *8, parameter :: logfactor_inv = 0.455119613313419D0
      real *8, parameter :: c16           = 0.16D0
      real *8, parameter :: half          = 0.5D0

      real *8 :: r_center(3)
      real *8 :: r0mr1, r0mr2, r0mr3
      real *8 :: pan_len_sq, pan_len_l_sq, pan_len_r_sq
      real *8 :: tplus, tminus

      if (t_root.ge.1.0D0) then
         r_center(1) = rr(1)
         r_center(2) = rr(2)
         r_center(3) = rr(3)
      else if (t_root.le.-1.0D0) then
         r_center(1) = rl(1)
         r_center(2) = rl(2)
         r_center(3) = rl(3)
      else
         r_center(1) = r_root(1)
         r_center(2) = r_root(2)
         r_center(3) = r_root(3)
      end if

      r0mr1 = r0(1) - r_center(1)
      r0mr2 = r0(2) - r_center(2)
      r0mr3 = r0(3) - r_center(3)

      sqn_dist2 = c16 * (r0mr1*r0mr1 + r0mr2*r0mr2 + r0mr3*r0mr3)

      pan_len_sq = pan_len * pan_len

      if (t_root.ge.1.0D0) then
         pan_len_l_sq = factor_sq * pan_len_sq
         lenl = floor( log(pan_len_l_sq / sqn_dist2) * logfactor_inv )+1
         if (lenl.lt.1) lenl = 1

         lenr = 0
         len  = lenl + 1

      else if (t_root.le.-1.0D0) then
         pan_len_r_sq = factor_sq * pan_len_sq
         lenr = floor( log(pan_len_r_sq / sqn_dist2) * logfactor_inv )+1
         if (lenr.lt.1) lenr = 1

         lenl = 0
         len  = lenr + 1

      else
         tplus  = (t_root + 1.0D0)
         tminus = (1.0D0 - t_root)

         pan_len_l_sq = (half*factor)**2 * tplus*tplus  * pan_len_sq
         pan_len_r_sq = (half*factor)**2 * tminus*tminus * pan_len_sq

         lenl = floor( log(pan_len_l_sq / sqn_dist2) * logfactor_inv )+1
         if (lenl.lt.1) lenl = 1

         lenr = floor( log(pan_len_r_sq / sqn_dist2) * logfactor_inv )+1
         if (lenr.lt.1) lenr = 1

         len  = lenl + lenr + 1
      end if

      end subroutine lenofpantendforsaf


      subroutine line3adaptivernearr0forsaf(t_root, p, tgl, wgl, r,
     1               w_bclag, len, lenl, lenr, rp, w, n,
     2               r_up, rp_up, Br, w_up)
      implicit none
      integer *8, intent(in) :: p, n, len, lenl, lenr
      real *8, intent(in) :: t_root, tgl(p), wgl(p), w_bclag(p), w(p)
      real *8, intent(in) :: r(3,p), rp(3,p)
      real *8, intent(inout) :: r_up(3,n), rp_up(3,n), Br(p,n), w_up(n)

      integer *8 :: k, j, i, npl
      real *8 :: tgl_up(n), pan_t_mid, pan_t_len
      real *8 :: B(n,p), winv(p)
      real *8 :: pan_t_end(len), pan_t_end1(lenl+1), pan_t_end2(lenr+1)
      real *8 :: factor_inv, xdiff, temp, factor
      real *8 :: denom(n)
      integer *8 :: idx0, nd
      real *8 :: rrp(p,6), rrp_up(n,6)
      real *8 :: wj, tj, wi
      real *8 :: nx, ny, nz

      nd        = 3
      factor    = 3.0d0
      factor_inv = 1.0d0 / factor

      pan_t_end1 = 0.0d0
      pan_t_end2 = 0.0d0

      if (t_root.ge.1.0d0) then
        pan_t_end1(lenl+1) = 1.0d0
        do k = lenl, 2, -1
          pan_t_end1(k) = factor_inv * pan_t_end1(k+1)
        end do
        pan_t_end = 1.0d0 - 2.0d0 * pan_t_end1(lenl+1:1:-1)
        npl = lenl + 1

      else if (t_root.le.-1.0d0) then
        pan_t_end2(lenr+1) = 1.0d0
        do k = lenr, 2, -1
          pan_t_end2(k) = factor_inv * pan_t_end2(k+1)
        end do
        pan_t_end = 2.0d0 * pan_t_end2 - 1.0d0
        npl = 0

      else
        pan_t_end1(lenl+1) = 1.0d0
        do k = lenl, 2, -1
          pan_t_end1(k) = factor_inv * pan_t_end1(k+1)
        end do
        pan_t_end2(lenr+1) = 1.0d0
        do k = lenr, 2, -1
          pan_t_end2(k) = factor_inv * pan_t_end2(k+1)
        end do
        do k = lenl+1, 1, -1
          pan_t_end(k) = -1.0d0 + (t_root+1.0d0) *
     1                  (1.0d0 - pan_t_end1(lenl+2-k))
        end do
        do k = 1, lenr+1
          pan_t_end(lenl+k) = 1.0d0 + (1.0d0-t_root) *
     1                        (-1.0d0 + pan_t_end2(k))
        end do
        npl = lenl + 1
      end if

      do j = 1, len-1
        pan_t_mid = 0.5d0 * (pan_t_end(j) + pan_t_end(j+1))
        pan_t_len = 0.5d0 * (pan_t_end(j+1) - pan_t_end(j))
        idx0 = (j-1)*p
        do k = 1, p
          tgl_up(idx0+k) = pan_t_mid + tgl(k)*pan_t_len
          w_up(idx0+k)   = wgl(k)   * pan_t_len
        end do
      end do

      do i = 1, n
        denom(i) = 0.0d0
      end do

      do j = 1, p
        wj = w_bclag(j)
        tj = tgl(j)
        do i = 1, n
          xdiff = tgl_up(i) - tj
          temp  = wj / xdiff
          B(i,j)   = temp
          denom(i) = denom(i) + temp
        end do
      end do
      do i = 1, n
        denom(i) = 1.0d0 / denom(i)
      end do
      do j = 1, p
        do i = 1, n
          B(i,j) = B(i,j) * denom(i)
        end do
      end do

      do k = 1, p
        rrp(k,1) = r(1,k)
        rrp(k,2) = r(2,k)
        rrp(k,3) = r(3,k)
        rrp(k,4) = rp(1,k)
        rrp(k,5) = rp(2,k)
        rrp(k,6) = rp(3,k)
      end do

      call dmatmat(n, p, B, 2*nd, rrp, rrp_up)

      do i = 1, n
        r_up(1,i)  = rrp_up(i,1)
        r_up(2,i)  = rrp_up(i,2)
        r_up(3,i)  = rrp_up(i,3)
        rp_up(1,i) = rrp_up(i,4)
        rp_up(2,i) = rrp_up(i,5)
        rp_up(3,i) = rrp_up(i,6)
      end do

      do i = 1, n
        nx = rp_up(1,i)
        ny = rp_up(2,i)
        nz = rp_up(3,i)
        w_up(i) = sqrt(nx*nx + ny*ny + nz*nz) * w_up(i)
      end do

      do k = 1, p
        winv(k) = 1.0d0 / w(k)
      end do

      do k = 1, p
        wi = winv(k)
        do i = 1, n
          B(i,k) = wi * w_up(i) * B(i,k)
        end do
      end do

      do j = 1, p
        do i = 1, n
          Br(j,i) = B(i,j)
        end do
      end do

      end subroutine line3adaptivernearr0forsaf

      subroutine line3nearrootf_real(tj, A, Legmat, wj, xj, yj, zj, n,
     1                         X, Y, Z, rho, rfc,
     2             troot_real,xroot_real,yroot_real,zroot_real,rat1temp)
      implicit none
      integer *8, intent(in) :: n
      real *8, intent(in) :: tj(n), wj(n), xj(n), yj(n), zj(n), rho
      real *8, intent(in) :: A(n,n), Legmat(n,n)
      real *8, intent(in) :: X, Y, Z
      real *8, intent(in) :: rat1temp(2,n)
      real *8, intent(inout) :: troot_real, xroot_real
      real *8, intent(inout) :: yroot_real, zroot_real
      integer *8, intent(inout) :: rfc

      complex *16 :: troot, xroot, yroot, zroot
      integer, parameter :: nmax_expa = 16
      integer *8 :: n_expa, i, j, converged
      real *8 :: xhat(nmax_expa), yhat(nmax_expa), zhat(nmax_expa)
      real *8 :: rat1(2, nmax_expa)
      real *8 :: P(nmax_expa)

      complex *16 :: tinit
      logical :: cp

      n_expa = min(nmax_expa, n)

      do i = 1, n_expa
         xhat(i) = 0.0d0
         yhat(i) = 0.0d0
         zhat(i) = 0.0d0
         do j = 1, n
            xhat(i) = xhat(i) + Legmat(i,j)*xj(j)
            yhat(i) = yhat(i) + Legmat(i,j)*yj(j)
            zhat(i) = zhat(i) + Legmat(i,j)*zj(j)
         end do
      end do

      do i = 1, n_expa
         rat1(1,i) = rat1temp(1,i)
         rat1(2,i) = rat1temp(2,i)
      end do

      tinit = CMPLX(0.0d0, 0.0d0, KIND=16)
      call rootfinderinitialguess(tj, xj, yj, zj, n, X, Y, Z,
     1                            tinit)

      cp = (abs(tinit + sqrt(tinit-1.0d0)*sqrt(tinit+1.0d0))
     1        < 1.5d0*rho)

      rfc = int(0, KIND=8)
      troot = CMPLX(0.0d0, 0.0d0, KIND=16)
      troot_real = 0.0d0
      xroot_real = 0.0d0
      yroot_real = 0.0d0
      zroot_real = 0.0d0

      if (cp) then
        converged = int(0, KIND=8)
        call rootfinderf_newton(xhat, yhat, zhat, n_expa, X, Y, Z,
     1                   tinit, troot, converged, rat1)

        troot_real = REAL(troot)   
        call legendrescalarf_real(n_expa-1, troot_real, P, rat1)

        ! 
        xroot_real = 0.0d0
        yroot_real = 0.0d0
        zroot_real = 0.0d0
        do i = 1, n_expa
           xroot_real = xroot_real + P(i)*xhat(i)
           yroot_real = yroot_real + P(i)*yhat(i)
           zroot_real = zroot_real + P(i)*zhat(i)
        end do

        if ( (converged == int(1, KIND=8)) .and.
     1   (abs(troot + sqrt(troot-1.0d0)*sqrt(troot+1.0d0)) < rho) ) then
          rfc = converged
        end if
      end if

      end subroutine line3nearrootf_real

      subroutine rootfinderf_newton(xhat,yhat,zhat,n,x0,y0,z0,tinit,
     1                      troot,converged,rat1)
      implicit none
      integer *8, intent(in) :: n
      real *8, intent(in) :: xhat(n), yhat(n), zhat(n), x0, y0, z0
      complex *16, intent(in) :: tinit
      complex *16, intent(inout) :: troot
      integer *8, intent(inout) :: converged
      real *8, intent(in) :: rat1(2,n)
      ! real *8, intent(inout) :: absres
      real *8 :: tol, absres
      integer *8 :: maxiter_newton, maxiter_muller, i, iter
      complex *16 :: F, Fprime, t, Fp, tp, Fpp, tpp, dt
      complex *16 :: P(n), D(n)
      complex *16 :: x, xp, y, yp, z, zp, dx, dy, dz, q, A, B, C, d1, d2

      t = tinit
      tol = 1.0d-14
      converged = 0
      maxiter_newton = 20
      maxiter_muller = 20

      ! Newton
      do iter = 1, maxiter_newton
        call legendrederivscalarf(n-1,t,P,D,rat1)
        x  = (0.0d0,0.0d0)
        y  = (0.0d0,0.0d0)
        z  = (0.0d0,0.0d0)
        xp = (0.0d0,0.0d0)
        yp = (0.0d0,0.0d0)
        zp = (0.0d0,0.0d0)
        do i = 1, n
          x  = x  + P(i) * xhat(i)
          y  = y  + P(i) * yhat(i)
          z  = z  + P(i) * zhat(i)
          xp = xp + D(i) * xhat(i)
          yp = yp + D(i) * yhat(i)
          zp = zp + D(i) * zhat(i)
        end do
        ! compute F and F'
        dx = x - x0
        dy = y - y0
        dz = z - z0
        F = dx*dx + dy*dy + dz*dz
        Fprime = 2*(dx*xp + dy*yp + dz*zp)
        dt = -F/Fprime
        ! update history
        tpp = tp
        Fpp = Fp
        Fp = F
        tp = t
        ! update root
        t = t + dt
        absres = abs(dt)
        if (absres < tol) then
          converged = 1
          exit
        endif
      enddo 
      
      if (converged == 1) then
        troot = t
        return
      endif

      end subroutine rootfinderf_newton


      subroutine legendrescalarf_real(n,x,P,rat1)
      implicit none
      integer *8, intent(in) :: n
      real *8, intent(in) :: x
      real *8, intent(inout), dimension(n+1) :: P
      real *8, intent(in), dimension(2,n+1) :: rat1
      integer *8 :: l
      P(1) = 1.0d0
      if (n >= 1) then
        P(2) = x
      endif
      do l=1,n-1
          P(l+2) = rat1(1,l+1)*x*P(l+1) + rat1(2,l+1)*P(l)
      enddo
      end subroutine legendrescalarf_real

c     omega_all fun again... 
      subroutine omegaall0123vrnewf(m, tx, m_all, nbd, nterms,    
     1                     morder, onm0, onm1, onm2, onm3, hdim, ijIdx, 
     2                     omega0, omega1, omega2, omega3)
c     
*/    
      implicit none
      integer *8, intent(in) :: m, nbd, hdim, morder, nterms
      real *8, intent(in), dimension(3,m) :: tx
      real *8, intent(in), dimension(morder,nbd,m) :: m_all
      real *8, intent(in), dimension(4,hdim,nbd) :: onm0, onm1 
      real *8, intent(in), dimension(4,hdim,nbd) :: onm2, onm3
      integer *8, intent(in), dimension(2,hdim) :: ijIdx
      
      integer *8, dimension(hdim) :: ijIdxsum
      integer *8 :: i, ii, j, k, l, ll, tmpidx, tmpidx2
      real *8, dimension(nbd*morder) :: mk_j
      real *8 :: pi, tmp
      real *8, intent(out), dimension(hdim,m) :: omega0,omega1
      real *8, intent(out), dimension(hdim,m) :: omega2,omega3
      real *8, dimension(4) :: mk_j_ii
      
      pi = 4.0d0*atan(1.0d0)
      tmp = -1.0d0/(4.0d0*pi)
      omega0 = 0.0D0
      omega1 = 0.0D0
      omega2 = 0.0D0
      omega3 = 0.0D0

c     loop over target      
      do j=1,m ! O(p^3) per target
c       integral along the edge, loop over edge quadrature node        
        do k=1,nbd ! O(3*p)
          i = 0
c         loop over moment nterms
          do ii=1,nterms ! O(p nterms)
            tmpidx = ii+1+(k-1)*morder
c           1 nterms higher moment            
            mk_j_ii(1) = tmp*m_all(ii+2,k,j)
c           some moment term depends on target linearly
            mk_j_ii(2) = tx(1,j)*tmp*m_all(ii+1,k,j)
            mk_j_ii(3) = tx(2,j)*tmp*m_all(ii+1,k,j)
            mk_j_ii(4) = tx(3,j)*tmp*m_all(ii+1,k,j)
c           loop over basis of that same moment nterms            
            do l=1,ii ! O(ii<=p)
c             basis count increases              
              i = i + 1 
              omega0(i,j) = omega0(i,j) 
     1                    + dot_product(mk_j_ii,onm0(:,i,k))
              omega1(i,j) = omega1(i,j) 
     1                    + dot_product(mk_j_ii,onm1(:,i,k))
              omega2(i,j) = omega2(i,j) 
     1                    + dot_product(mk_j_ii,onm2(:,i,k))
              omega3(i,j) = omega3(i,j) 
     1                    + dot_product(mk_j_ii,onm3(:,i,k))
            enddo
          enddo
        enddo
      enddo 

      end subroutine omegaall0123vrnewf  

      
