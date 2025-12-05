      ! -----------------------------------------------------------------
      subroutine Lap3dSLP_closepaneladp_vr_guru(m,tx,
     1                                     n,sx,sw,snx,iside,Aslp,Adlp)
      implicit none
      integer *8, intent(in) :: m, n, iside
      real *8, intent(in) :: tx(3,m), sx(3,n), sw(n), snx(3,n)
      real *8, intent(inout) :: Aslp(m,n), Adlp(m,n)
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
      call Lap3dSLP_closepaneladp_vr(m,tx,nterms,
     1                      n,sx,sw,snx,umatr,iside,hdim,
     2                      alpha,sxc,r_vert,
     3                      sbdnp,len,nbd,nquad,sxbd,stangbd,swbd,sspbd,
     4                      tgl,wgl,w_bclag,Dgl,Agl,Legmat,bclagmatlr,
     5                      onm0,onm1,onm2,onm3,Fbd,Fxbd,Fybd,Fzbd,
     6                      Mmatrix,isimd,Aslp,Adlp,Omega)

      deallocate(sxbd, stangbd)
      deallocate(swbd, sspbd)
      deallocate(tgl, wgl, w_bclag)
      deallocate(Dgl, Agl)
      deallocate(Legmat, bclagmatlr)
      deallocate(onm0, onm1, onm2, onm3)
      deallocate(Fbd, Fxbd, Fybd, Fzbd)
      end subroutine Lap3dSLP_closepaneladp_vr_guru

c     ... 
      subroutine Lap3dSLP_closepaneladp_vr(m,tx,nterms,
     1                      n,sx,sw,snx,umatr,iside,hdim,
     2                      alpha,sxc,r_vert,
     3                      sbdnp,len,nbd,nquad,sxbd,stangbd,swbd,sspbd,
     4                      tgl,wgl,w_bclag,Dgl,Agl,Legmat,bclagmatlr,
     5                      onm0,onm1,onm2,onm3,Fbd,Fxbd,Fybd,Fzbd,
     6                      Mmatrix,isimd,Aslp,Adlp,Omega)
      implicit none
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
      real *8, intent(inout) :: Mmatrix(4*n,4*n), Omega(4*n,m)
      real *8, intent(inout) :: Adlp(m,n), Aslp(m,n)

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
      real *8 :: Fbd_real(nbd,n), F1bd(nbd,n), F2bd(nbd,n), F3bd(nbd,n)
      integer *8 :: morder
      real *8, allocatable :: m_all_adp(:,:,:), n_all_adp(:,:,:)
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
      real *8 :: mk_j(nquad,nterms+3), nk_j(nquad,nterms+3)
      integer *8 :: ijIdx(2,n), idx, idx_k_end
      real *8 :: omega0slp(n,m),omega1slp(n,m),omega2slp(n,m),
     1           omega3slp(n,m)
      real *8 :: omega0dlp(n,m),omega1dlp(n,m),omega2dlp(n,m),
     1           omega3dlp(n,m)
      complex *16 :: Fc(n,(nterms+1)*(nterms+1))
      complex *16 :: Fx_c(n,(nterms+1)*(nterms+1))
      complex *16 :: Fy_c(n,(nterms+1)*(nterms+1))
      complex *16 :: Fz_c(n,(nterms+1)*(nterms+1))
      integer *8 :: info
      real *8 :: F0sx(n,n), F1sx(n,n), F2sx(n,n), F3sx(n,n)
      real *8 :: MmatrixT(4*n,4*n), Atmp(4*n,m)
      real *8 :: onm0slp(5,n,nbd), onm1slp(5,n,nbd)
      real *8 :: onm2slp(5,n,nbd), onm3slp(5,n,nbd)
      real *8 :: Omega_slp(m,4*n)
      real *8 :: Amu_c(4*n,n), rhs_basis(4*n,n), dlm_basis(4*n,n)
      real *8 :: fval_closed_basis(m,n), fval_added_basis(m,n)
      real *8 :: fval_basis(m,n)
      real *8 :: dlm_0_basis(n,n), dlm_1_basis(n,n)
      real *8 :: dlm_2_basis(n,n), dlm_3_basis(n,n)
      real *8 :: rho_basis(4*n,n), glm_basis(4*n,n)

      pi = 4.0d0*atan(1.0d0)
      o_offset = 0.063d0
      hdim = nterms*(nterms+1)/2

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

      korder = nterms-1
      vmatr = 0.0d0
      umatr = 0.0d0
      call koorn_vals2coefs_coefs2vals(korder,n,umatr,vmatr)

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

      Omega = 0.0d0
      Adlp = 0.0d0
      Aslp = 0.0d0
      Fbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fxbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fybd = CMPLX((0.0d0, 0.0d0), KIND=16)
      Fzbd = CMPLX((0.0d0, 0.0d0), KIND=16)
      call l3dtavecevalmatf(sxbd,nbd,nterms,Fbd,Fxbd,Fybd,Fzbd,ier)
      do k = 1,n
        Fbd_real(:,k) = dble(Fbd(:,idxvec(k)))
        F1bd(:,k) = dble(Fxbd(:,idxvec(k)))
        F2bd(:,k) = dble(Fybd(:,idxvec(k)))
        F3bd(:,k) = dble(Fzbd(:,idxvec(k)))
      enddo

      onm0slp = 0.0d0
      onm1slp = 0.0d0
      onm2slp = 0.0d0
      onm3slp = 0.0d0
      onm0 = 0.0d0
      onm1 = 0.0d0
      onm2 = 0.0d0
      onm3 = 0.0d0
      call omeganmslp_precomp(nbd,hdim,Fbd_real,F1bd,F2bd,F3bd,
     1     sxbd,stangbd,swbd,onm0slp,onm1slp,onm2slp,onm3slp,
     2     onm0,onm1,onm2,onm3)

      morder = nterms+3
      allocate( m_all_adp(morder,nbd,m))
      allocate( n_all_adp(morder,nbd,m))
      m_all_adp = 0.0d0
      n_all_adp = 0.0d0

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
        do j = 1,m
          r0j = txnew(:, j)
          x0 = r0j(1)
          y0 = r0j(2)
          z0 = r0j(3)
          rfc = int(0, KIND=8)
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
            call lenofpantendforsaf(nquad,tgl,wgl,r_ell,rr,rl,r0j,
     1              troot_real,r_root,w_bclag,Dgl,sqn_dist,
     2              pan_len,lenj,lenl,lenr,rp,w)
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
            flag = int(0, KIND=8)
            allocate(mk_j_up(nquad_up,morder))
            allocate(nk_j_up(nquad_up,morder))
            mk_j_up = 0.0d0
            nk_j_up = 0.0d0
            call momentsad_vr(r0j, nquad_up, r_up, nterms+2, flag, 
     1              nk_j_up, mk_j_up)
            mk_j = matmul(Br, mk_j_up)
            nk_j = matmul(Br, nk_j_up)

            deallocate(r_up,rp_up)
            deallocate(Br,w_up)
            deallocate(mk_j_up)
            deallocate(nk_j_up)
          else
            flag = int(0, KIND=8)
            mk_j = 0.0d0 
            nk_j = 0.0d0
            call momentsad_vr(r0j, nquad, r_ell, nterms+2, flag, 
     1              nk_j, mk_j)
          end if  
          m_all_adp(1:morder,idx_ell_start:idx_ell_end,j) = 
     1              transpose(mk_j)
          n_all_adp(1:morder,idx_ell_start:idx_ell_end,j) = 
     1              transpose(nk_j)
        enddo
      enddo

      ijIdx = int(0, KIND=8)
      idx = int(1, KIND=8)
      do k=1,nterms
        idx_k_end = idx + k - 1
        ijIdx(1, idx:idx_k_end) = ijIdx(1, idx:idx_k_end) + (k - 1)
        idx = idx + k
      enddo
      omega0slp = 0.0d0
      omega1slp = 0.0d0
      omega2slp = 0.0d0
      omega3slp = 0.0d0
      omega0dlp = 0.0d0
      omega1dlp = 0.0d0
      omega2dlp = 0.0d0
      omega3dlp = 0.0d0
      call omegasdlpall0123vrnewf(m, txnew, n_all_adp, m_all_adp, nbd, 
     1    nterms, morder, onm0slp, onm1slp, onm2slp, onm3slp, 
     2    onm0, onm1, onm2, onm3, hdim, ijIdx, omega0slp, omega1slp, 
     3    omega2slp, omega3slp, omega0dlp, omega1dlp, omega2dlp, 
     4    omega3dlp)

      Omega_slp = 0.0d0
      do k=1,n
         Omega_slp(:,k) = -omega0slp(k,:)/alpha
         Omega_slp(:,n+k) = omega1slp(k,:)/alpha
         Omega_slp(:,2*n+k) = omega2slp(k,:)/alpha
         Omega_slp(:,3*n+k) = omega3slp(k,:)/alpha
      enddo

      do k=1,m
         Omega(1:n,k) = omega0dlp(1:n,k)
         Omega(n+1:2*n,k) = -omega1dlp(1:n,k)
         Omega(2*n+1:3*n,k) = -omega2dlp(1:n,k)
         Omega(3*n+1:4*n,k) = -omega3dlp(1:n,k)
      end do

      Fc = (0.0d0,0.0d0)
      Fx_c = (0.0d0,0.0d0)
      Fy_c = (0.0d0,0.0d0)
      Fz_c = (0.0d0,0.0d0)
      info = 0
      call l3dtavecevalmatf(sxnew, n, nterms, Fc, Fx_c, Fy_c, Fz_c,info)

      F0sx = 0.0d0
      do k=1,hdim
        F0sx(:,k) = dble(Fc(:,idxvec(k)))
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
      Adlp = transpose(Atmp(1:n,:))

      Amu_c = 0.0d0
      do k=1,n
         Amu_c(n+k,k) = snxnew(1,k)
         Amu_c(2*n+k,k) = snxnew(2,k)
         Amu_c(3*n+k,k) = snxnew(3,k)
      enddo

      dlm_basis = matmul(Amu_c,vmatr)
      call dgausselimvec(4*n,Mmatrix,n,dlm_basis,info)

      fval_closed_basis = matmul(Omega_slp,dlm_basis)
      dlm_0_basis = (dlm_basis(1:n,:))
      dlm_1_basis = (dlm_basis(n+1:2*n,:))
      dlm_2_basis = (dlm_basis(2*n+1:3*n,:))
      dlm_3_basis = (dlm_basis(3*n+1:4*n,:))

      rho_basis = 0.0d0
      rho_basis(1:n,:) = -matmul(F0sx,dlm_0_basis)
      rho_basis(n+1:2*n,:) = -matmul(F0sx,dlm_1_basis)
      rho_basis(2*n+1:3*n,:) = -matmul(F0sx,dlm_2_basis)
      rho_basis(3*n+1:4*n,:) = -matmul(F0sx,dlm_3_basis)

      glm_basis = rho_basis
      call dgausselimvec(4*n,Mmatrix,n,glm_basis,info)
      fval_added_basis = matmul(transpose(Omega),glm_basis)
      fval_basis = fval_closed_basis - fval_added_basis/alpha
      Aslp = matmul(fval_basis,umatr)

      deallocate( m_all_adp )
      deallocate( n_all_adp )

      end subroutine Lap3dSLP_closepaneladp_vr

c     omega_all fun again... 
      subroutine omegasdlpall0123vrnewf(m, tx, n_all, m_all, nbd,nterms,    
     1               morder, onm0slp, onm1slp, onm2slp, onm3slp, 
     2               onm0dlp, onm1dlp, onm2dlp, onm3dlp, hdim, ijIdx, 
     3               omega0slp, omega1slp, omega2slp, omega3slp, 
     4               omega0dlp, omega1dlp, omega2dlp, omega3dlp)
c     
*/    
      implicit none
      integer *8, intent(in) :: m, nbd, hdim, morder, nterms
      real *8, intent(in), dimension(3,m) :: tx
      real *8, intent(in), dimension(morder,nbd,m) :: m_all, n_all
      real *8, intent(in), dimension(5,hdim,nbd) :: onm0slp, onm1slp 
      real *8, intent(in), dimension(5,hdim,nbd) :: onm2slp, onm3slp
      real *8, intent(in), dimension(4,hdim,nbd) :: onm0dlp, onm1dlp 
      real *8, intent(in), dimension(4,hdim,nbd) :: onm2dlp, onm3dlp
      integer *8, intent(in), dimension(2,hdim) :: ijIdx
      
      integer *8, dimension(hdim) :: ijIdxsum
      integer *8 :: i, ii, j, k, l, ll, tmpidx, tmpidx2
      real *8, dimension(nbd*morder) :: mk_j
      real *8 :: pi, tmp
      real *8, intent(out), dimension(hdim,m) :: omega0slp,omega1slp
      real *8, intent(out), dimension(hdim,m) :: omega2slp,omega3slp
      real *8, intent(out), dimension(hdim,m) :: omega0dlp,omega1dlp
      real *8, intent(out), dimension(hdim,m) :: omega2dlp,omega3dlp
      real *8, dimension(4) :: mk_j_ii_dlp
      real *8, dimension(5) :: mk_j_ii_slp
      
      pi = 4.0d0*atan(1.0d0)
      tmp = -1.0d0/(4.0d0*pi)
      omega0slp = 0.0D0
      omega1slp = 0.0D0
      omega2slp = 0.0D0
      omega3slp = 0.0D0
      omega0dlp = 0.0D0
      omega1dlp = 0.0D0
      omega2dlp = 0.0D0
      omega3dlp = 0.0D0

c     loop over target      
      do j=1,m ! O(p^3) per target
c       integral along the edge, loop over edge quadrature node        
        do k=1,nbd ! O(3*p)
          i = 0
c         loop over moment nterms
          do ii=1,nterms ! O(p nterms)
c           1 nterms higher moment            
            mk_j_ii_dlp(1) = tmp*m_all(ii+2,k,j)
c           some moment term depends on target linearly
            mk_j_ii_dlp(2) = tx(1,j)*tmp*m_all(ii+1,k,j)
            mk_j_ii_dlp(3) = tx(2,j)*tmp*m_all(ii+1,k,j)
            mk_j_ii_dlp(4) = tx(3,j)*tmp*m_all(ii+1,k,j)
            ! 
            mk_j_ii_slp(1) = tmp*m_all(ii+3,k,j)
            ! 
            mk_j_ii_slp(2) = tx(1,j)*tmp*m_all(ii+2,k,j)
            mk_j_ii_slp(3) = tx(2,j)*tmp*m_all(ii+2,k,j)
            mk_j_ii_slp(4) = tx(3,j)*tmp*m_all(ii+2,k,j)
            ! 
            mk_j_ii_slp(5) = tmp*n_all(ii+1,k,j)
c           loop over basis of that same moment nterms            
            do l=1,ii ! O(ii<=p)
c             basis count increases              
              i = i + 1 
            !   dlp
              omega0dlp(i,j) = omega0dlp(i,j) 
     1                    + dot_product(mk_j_ii_dlp,onm0dlp(:,i,k))
              omega1dlp(i,j) = omega1dlp(i,j) 
     1                    + dot_product(mk_j_ii_dlp,onm1dlp(:,i,k))
              omega2dlp(i,j) = omega2dlp(i,j) 
     1                    + dot_product(mk_j_ii_dlp,onm2dlp(:,i,k))
              omega3dlp(i,j) = omega3dlp(i,j) 
     1                    + dot_product(mk_j_ii_dlp,onm3dlp(:,i,k))
            !   slp
              omega0slp(i,j) = omega0slp(i,j) 
     1                    + dot_product(mk_j_ii_slp,onm0slp(:,i,k))
              omega1slp(i,j) = omega1slp(i,j) 
     1                    + dot_product(mk_j_ii_slp,onm1slp(:,i,k))
              omega2slp(i,j) = omega2slp(i,j) 
     1                    + dot_product(mk_j_ii_slp,onm2slp(:,i,k))
              omega3slp(i,j) = omega3slp(i,j) 
     1                    + dot_product(mk_j_ii_slp,onm3slp(:,i,k))
            enddo
          enddo
        enddo
      enddo 

      end subroutine omegasdlpall0123vrnewf


c     Precompute onm tensors (fast coef version)
      subroutine omeganmslp_precompff(nbd, hdim, F, F1, F2, F3,
     1     sxbd, stangbd, swbd, onm0slp, onm1slp, onm2slp, onm3slp,
     2     onm0dlp, onm1dlp, onm2dlp, onm3dlp)
      implicit none
      integer *8, intent(in) :: nbd, hdim
      real *8, intent(in) :: F(nbd,hdim)
      real *8, intent(in) :: F1(nbd,hdim), F2(nbd,hdim), F3(nbd,hdim)
      real *8, intent(in) :: sxbd(3,nbd), stangbd(3,nbd), swbd(nbd)
      real *8, intent(inout) :: onm0slp(5,hdim,nbd), onm1slp(5,hdim,nbd)
      real *8, intent(inout) :: onm2slp(5,hdim,nbd), onm3slp(5,hdim,nbd)
      real *8, intent(inout) :: onm0dlp(4,hdim,nbd), onm1dlp(4,hdim,nbd)
      real *8, intent(inout) :: onm2dlp(4,hdim,nbd), onm3dlp(4,hdim,nbd)

      integer *8 :: i, j, k
      real *8 :: dx(nbd), dy(nbd), dz(nbd)
      real *8 :: sx1, sx2, sx3, tx, ty, tz
      real *8 :: ft, f1t, f2t, f3t, mr_dot_f
      real *8 :: c1, c2, c3, c4, c5, c6
      real *8 :: acoef, bcoef, ccoef
      real *8 :: q01d(4), q02d(4), q03d(4)
      real *8 :: q11d(4), q12d(4), q13d(4)
      real *8 :: q21d(4), q22d(4), q23d(4)
      real *8 :: q31d(4), q32d(4), q33d(4)
      real *8 :: q01s(5), q02s(5), q03s(5)
      real *8 :: q11s(5), q12s(5), q13s(5)
      real *8 :: q21s(5), q22s(5), q23s(5)
      real *8 :: q31s(5), q32s(5), q33s(5)

      do k = 1, nbd
         dx(k) = stangbd(1,k)*swbd(k)
         dy(k) = stangbd(2,k)*swbd(k)
         dz(k) = stangbd(3,k)*swbd(k)
      end do

      onm0slp = 0.0d0
      onm1slp = 0.0d0
      onm2slp = 0.0d0
      onm3slp = 0.0d0
      onm0dlp = 0.0d0
      onm1dlp = 0.0d0
      onm2dlp = 0.0d0
      onm3dlp = 0.0d0

      do k = 1, nbd
         sx1 = sxbd(1,k)
         sx2 = sxbd(2,k)
         sx3 = sxbd(3,k)
         tx  = dx(k)
         ty  = dy(k)
         tz  = dz(k)

         c1 = tx*sx3
         c2 = tx*sx2
         c3 = ty*sx1
         c4 = ty*sx3
         c5 = tz*sx2
         c6 = tz*sx1

         acoef = -c4 + c5
         bcoef =  c1 - c6
         ccoef = -c2 + c3

         do i = 1, hdim
            ft  = F(k,i)
            f1t = F1(k,i)
            f2t = F2(k,i)
            f3t = F3(k,i)

            mr_dot_f = -(sx1*f1t + sx2*f2t + sx3*f3t)

            ! ----- DLP qnm -----
            q01d(1) = -sx2*f3t + sx3*f2t
            q01d(2) =  0.0d0
            q01d(3) =  f3t
            q01d(4) = -f2t

            q02d(1) = -sx3*f1t + sx1*f3t
            q02d(2) = -f3t
            q02d(3) =  0.0d0
            q02d(4) =  f1t

            q03d(1) = -sx1*f2t + sx2*f1t
            q03d(2) =  f2t
            q03d(3) = -f1t
            q03d(4) =  0.0d0

            q11d(1) =  mr_dot_f + 2.0d0*sx1*f1t
            q11d(2) = -f1t
            q11d(3) =  f2t
            q11d(4) =  f3t

            q12d(1) =  sx1*f2t + sx2*f1t
            q12d(2) = -f2t
            q12d(3) = -f1t
            q12d(4) =  0.0d0

            q13d(1) =  sx1*f3t + sx3*f1t
            q13d(2) = -f3t
            q13d(3) =  0.0d0
            q13d(4) = -f1t

            q21d(1) =  sx2*f1t + sx1*f2t
            q21d(2) = -f2t
            q21d(3) = -f1t
            q21d(4) =  0.0d0

            q22d(1) =  mr_dot_f + 2.0d0*sx2*f2t
            q22d(2) =  f1t
            q22d(3) = -f2t
            q22d(4) =  f3t

            q23d(1) =  sx2*f3t + sx3*f2t
            q23d(2) =  0.0d0
            q23d(3) = -f3t
            q23d(4) = -f2t

            q31d(1) =  sx3*f1t + sx1*f3t
            q31d(2) = -f3t
            q31d(3) =  0.0d0
            q31d(4) = -f1t

            q32d(1) =  sx3*f2t + sx2*f3t
            q32d(2) =  0.0d0
            q32d(3) = -f3t
            q32d(4) = -f2t

            q33d(1) =  mr_dot_f + 2.0d0*sx3*f3t
            q33d(2) =  f1t
            q33d(3) =  f2t
            q33d(4) = -f3t

            ! ----- SLP qnm -----
            q01s = 0.0d0
            q02s = 0.0d0
            q03s = 0.0d0
            q11s = 0.0d0
            q12s = 0.0d0
            q13s = 0.0d0
            q21s = 0.0d0
            q22s = 0.0d0
            q23s = 0.0d0
            q31s = 0.0d0
            q32s = 0.0d0
            q33s = 0.0d0

            q01s(1) =  sx1*ft
            q01s(2) = -ft
            q01s(5) =  f1t

            q02s(1) =  sx2*ft
            q02s(3) = -ft
            q02s(5) =  f2t

            q03s(1) =  sx3*ft
            q03s(4) = -ft
            q03s(5) =  f3t

            q12s(1) =  sx3*ft
            q12s(4) = -ft
            q12s(5) = -f3t

            q13s(1) = -sx2*ft
            q13s(3) =  ft
            q13s(5) =  f2t

            q21s(1) = -sx3*ft
            q21s(4) =  ft
            q21s(5) =  f3t

            q23s(1) =  sx1*ft
            q23s(2) = -ft
            q23s(5) = -f1t

            q31s(1) =  sx2*ft
            q31s(3) = -ft
            q31s(5) = -f2t

            q32s(1) = -sx1*ft
            q32s(2) =  ft
            q32s(5) =  f1t

            do j = 1, 4
               onm0dlp(j,i,k) = acoef*q01d(j) + bcoef*q02d(j)
     1             + ccoef*q03d(j)
               onm1dlp(j,i,k) = acoef*q11d(j) + bcoef*q12d(j)
     1             + ccoef*q13d(j)
               onm2dlp(j,i,k) = acoef*q21d(j) + bcoef*q22d(j)
     1             + ccoef*q23d(j)
               onm3dlp(j,i,k) = acoef*q31d(j) + bcoef*q32d(j)
     1             + ccoef*q33d(j)

               onm0slp(j,i,k) = acoef*q01s(j) + bcoef*q02s(j)
     1             + ccoef*q03s(j)
               onm1slp(j,i,k) = acoef*q11s(j) + bcoef*q12s(j)
     1             + ccoef*q13s(j)
               onm2slp(j,i,k) = acoef*q21s(j) + bcoef*q22s(j)
     1             + ccoef*q23s(j)
               onm3slp(j,i,k) = acoef*q31s(j) + bcoef*q32s(j)
     1             + ccoef*q33s(j)
            end do

            j = 5
            onm0slp(j,i,k) = acoef*q01s(j) + bcoef*q02s(j)
     1          + ccoef*q03s(j)
            onm1slp(j,i,k) = acoef*q11s(j) + bcoef*q12s(j)
     1          + ccoef*q13s(j)
            onm2slp(j,i,k) = acoef*q21s(j) + bcoef*q22s(j)
     1          + ccoef*q23s(j)
            onm3slp(j,i,k) = acoef*q31s(j) + bcoef*q32s(j)
     1          + ccoef*q33s(j)
         end do
      end do

      end subroutine omeganmslp_precompff
