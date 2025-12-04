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

      
c     Precompute onm tensors for SLP/DLP
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
