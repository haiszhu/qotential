c     Fortran evalTensorProductHarmonicGrad fun. Hai 01/01/23
c     
      subroutine evaltensorproductharmonicgrad(nt,r,order,
     1           fx,fy,fz,f,ijidx)
c     
*/ 
      implicit none
      integer, intent(in) :: nt, order
      real *8, intent(in), dimension(3,nt) :: r
      real *8, intent(out), dimension(nt,order**2) :: fx, fy, fz, f
      integer, intent(out), dimension(2,order**2) :: ijidx
      
      integer :: k
      real *8, dimension(nt,order+1) :: rxpow, rypow
      real *8, dimension(nt,2*order+1) :: rzpow
      real *8, dimension(order**2,order**2) :: coeffs_all
      integer, dimension(order**2,order**2) ::powx_all,powy_all,powz_all
      integer, dimension(order**2) :: lens
      real *8, dimension(order**2) :: coeffs_ext
      integer, dimension(order**2) :: powx_ext, powy_ext, powz_ext
      integer :: lens_k
      real *8, dimension(nt) :: fx_k, fy_k, fz_k, f_k


      rxpow = 1.0d0
      rypow = 1.0d0
      rzpow = 1.0d0
      do k=1,order
        rxpow(:,k+1) = rxpow(:,k)*r(1,:)
        rypow(:,k+1) = rypow(:,k)*r(2,:)
      enddo
      do k=1,2*order
        rzpow(:,k+1) = rzpow(:,k)*r(3,:)
      enddo

      call hijcoeffsall(order,coeffs_all,powx_all,powy_all,powz_all,
     1                  lens,ijidx)

      fx = 0.0d0
      fy = 0.0d0
      fz = 0.0d0
      f = 0.0d0
! C$OMP PARALLEL DO DEFAULT(SHARED)
! C$OMP$PRIVATE(k,coeffs_ext,powx_ext,powy_ext,powz_ext,lens_k)
! C$OMP$PRIVATE(fx_k,fy_k,fz_k,f_k)
      do k=1,order**2
        coeffs_ext = coeffs_all(k,:)
        powx_ext = powx_all(k,:)
        powy_ext = powy_all(k,:)
        powz_ext = powz_all(k,:)
        lens_k = lens(k)

        call fxfyfzf(order**2, coeffs_ext, powx_ext, powy_ext, powz_ext,
     1           lens_k, nt, order, rxpow, rypow, rzpow, 
     2           fx_k, fy_k, fz_k, f_k)
        
        fx(:,k) = fx_k
        fy(:,k) = fy_k
        fz(:,k) = fz_k
        f(:,k) = f_k

      enddo
! C$OMP END PARALLEL DO
 
      end subroutine evaltensorproductharmonicgrad

c     Fortran fx(:,k) fy(:,k) fz(:,k) f(:,k) fun. Hai 01/01/23
c
      subroutine fxfyfzf(order2, coeffs, powx, powy, powz, lens_k,
     1           nt, order, rxpow, rypow, rzpow, 
     2           fx_k, fy_k, fz_k, f_k)
c     
*/ 
      implicit none
      integer, intent(in) :: order2, lens_k, nt, order
      real *8, intent(in), dimension(order2) :: coeffs
      integer, intent(in), dimension(order2) :: powx, powy, powz
      real *8, intent(in), dimension(nt,order+1) :: rxpow, rypow
      real *8, intent(in), dimension(nt,2*order+1) :: rzpow
      real *8, intent(out), dimension(nt) :: fx_k, fy_k, fz_k, f_k

      integer, dimension(3,lens_k) :: pow_k
      real *8, dimension(lens_k) :: coeffs_k
      integer, dimension(3,lens_k) :: powgradx, powgrady, powgradz
      integer :: j, l
      real *8, dimension(lens_k) :: coeffsgradx,coeffsgrady,coeffsgradz
      

      pow_k(1,:) = powx(1:lens_k)
      pow_k(2,:) = powy(1:lens_k)
      pow_k(3,:) = powz(1:lens_k)
      coeffs_k = coeffs(1:lens_k)

      powgradx = pow_k
      powgrady = pow_k
      powgradz = pow_k
      powgradx(1,:) = powgradx(1,:)-1
      powgrady(2,:) = powgrady(2,:)-1
      powgradz(3,:) = powgradz(3,:)-1
      do j=1,lens_k
        if (powgradx(1,j).eq.-1) then
          powgradx(1,j) = 0
        endif
        if (powgrady(2,j).eq.-1) then
          powgrady(2,j) = 0
        endif
        if (powgradz(3,j).eq.-1) then
          powgradz(3,j) = 0
        endif
      enddo

      coeffsgradx = coeffs_k*pow_k(1,:)
      coeffsgrady = coeffs_k*pow_k(2,:)
      coeffsgradz = coeffs_k*pow_k(3,:)

      fx_k = 0.0d0 
      fy_k = 0.0d0 
      fz_k = 0.0d0 
      f_k  = 0.0d0
      do j=1,lens_k
        fx_k = fx_k + coeffsgradx(j)*rxpow(:,powgradx(1,j)+1)
     1               *rypow(:,powgradx(2,j)+1)*rzpow(:,powgradx(3,j)+1)
        fy_k = fy_k + coeffsgrady(j)*rxpow(:,powgrady(1,j)+1)
     1               *rypow(:,powgrady(2,j)+1)*rzpow(:,powgrady(3,j)+1)
        fz_k = fz_k + coeffsgradz(j)*rxpow(:,powgradz(1,j)+1)
     1               *rypow(:,powgradz(2,j)+1)*rzpow(:,powgradz(3,j)+1)
        f_k = f_k + coeffs_k(j)*rxpow(:,pow_k(1,j)+1)
     1               *rypow(:,pow_k(2,j)+1)*rzpow(:,pow_k(3,j)+1)
      enddo


      end subroutine fxfyfzf

c     Fortran Hij_coeffs_all fun. Hai 01/01/23
      subroutine hijcoeffsall(order,coeffs,powx,powy,powz,lens,ijidx)
c     
*/ 
      implicit none
      integer, intent(in) :: order
      real *8, intent(out), dimension(order**2,order**2) :: coeffs
      integer, intent(out), dimension(order**2,order**2)::powx,powy,powz
      integer, intent(out),dimension(order**2) :: lens
      integer, intent(out), dimension(2,order**2) :: ijidx

      integer :: orderx, ordery, k, i, j, lens_k
      real *8, dimension(order**2) :: coeffs_k
      integer, dimension(3,order**2) :: pow_k

      orderx = order
      ordery = order

      call tdordering(orderx,ordery,ijidx)
      coeffs = 0.0d0
      powx = 0
      powy = 0
      powz = 0
      lens = 0
! C$OMP PARALLEL DO DEFAULT(SHARED)
! C$OMP$PRIVATE(k,i,j,lens_k,coeffs_k,pow_k)
      do k=1,order**2
        i = ijidx(1,k)
        j = ijidx(2,k)
        call hijcoeffs0(i,j,order**2,coeffs_k,pow_k,lens_k)
        lens(k) = lens_k
        coeffs(k,:) = coeffs_k
        powx(k,:) = pow_k(1,:)
        powy(k,:) = pow_k(2,:)
        powz(k,:) = pow_k(3,:)
      enddo
! C$OMP END PARALLEL DO 

      end subroutine hijcoeffsall

c     Fortran Hij_coeffs fun. Hai 01/01/23
c
      subroutine hijcoeffs0(i,j,orderxy,coeffs_ext,pow_ext,len)
c     
*/       
      implicit none
      integer, intent(in) :: i, j, orderxy
      real *8, intent(out), dimension(orderxy) :: coeffs_ext
      integer, intent(out), dimension(3,orderxy) :: pow_ext
      integer, intent(out) :: len

      real *8, dimension((i+1)*(j+1)) :: coeffs
      integer, dimension(3,(i+1)*(j+1)) :: pow
      integer :: l, k
      

      call hijcoeffs(i,j,coeffs,pow)
      len = size(coeffs)
      coeffs_ext = 0.0d0
      pow_ext = 0
      coeffs_ext((/(l,l=1,len)/)) = coeffs
      pow_ext(1,(/(l,l=1,len)/)) = pow(1,:)
      pow_ext(2,(/(l,l=1,len)/)) = pow(2,:)
      pow_ext(3,(/(l,l=1,len)/)) = pow(3,:)

      end subroutine hijcoeffs0

      

c     Fortran td_ordering fun. Hai 01/01/23
      subroutine tdordering(orderx,ordery,ijidx)
c     
*/ 
      implicit none
      integer, intent(in) :: orderx, ordery
      integer, intent(out), dimension(2,orderx*ordery) :: ijidx
      
      integer, dimension(orderx,ordery) :: iidx, jidx
      integer :: l, j, kstart, kend
      logical, dimension(orderx*ordery) :: tmpidx
      integer, dimension(orderx*ordery) :: iidx_rs, jidx_rs
      
      call meshgrid(ordery,(/(l,l=0,ordery-1)/),
     1              orderx,(/(l,l=0,orderx-1)/),jidx,iidx)
      ijIdx = 0
      iidx_rs = reshape(iidx, [orderx*ordery])
      jidx_rs = reshape(jidx, [orderx*ordery])

      kstart = 1
      do l=0,ordery+orderx-2
        tmpidx = (jidx_rs+iidx_rs).eq.l
        kend = kstart + count(tmpidx) - 1
        do j=1,orderx*ordery
          if (tmpidx(j)) then
            ijidx(1,kstart) = iidx_rs(j)
            ijidx(2,kstart) = jidx_rs(j)
            kstart = kstart + 1
          endif
        enddo
      enddo

      end subroutine tdordering

c     Fortran Hij_coeffs fun. Hai 01/01/23
c
      subroutine hijcoeffs(i,j,coeffs,pow)
c     
*/       
      implicit none
      integer, intent(in) :: i, j
      real *8, intent(out), dimension((i+1)*(j+1)) :: coeffs
      integer, intent(out), dimension(3,(i+1)*(j+1)) :: pow

      integer, dimension(i+1,j+1) :: xpow, ypow
      integer :: l, k, len
      real *8, dimension((i+1)*(j+1)) :: coeffs_updt
      

      call meshgrid(j+1,(/(l,l=0,j)/),i+1,(/(l,l=0,i)/),ypow,xpow)
      pow(1,:) = reshape(xpow, [(i+1)*(j+1)])
      pow(2,:) = reshape(ypow, [(i+1)*(j+1)])
      pow(3,:) = i+j+1 - (pow(1,:)+pow(2,:))

      coeffs = 0.0
      do k=0,floor(dble(i+j)/2)
        call lenofnckcoeffs0(i, j, k, len)
        call hijcoeffsk0(i, j, k, pow, len, coeffs_updt)
        coeffs = coeffs + coeffs_updt
      enddo


      end subroutine hijcoeffs

      subroutine meshgrid(n, x, m, y, xx, yy)
c     
      integer, intent(in) :: n, m
      integer, intent(in), dimension(n) :: x 
      integer, intent(in), dimension(m) :: y
      integer, intent(out), dimension(m,n) :: xx, yy
      
      xx = spread(x, 1, size(y))
      yy = spread(y, 2, size(x))

      end subroutine

c     Fortran fun for the update part of Hij_coeffs from \Delta^k
c     Hai 01/01/23
c
      subroutine hijcoeffsk0(i, j, k, pow, len, coeffs_updt)
c     
*/
      implicit none
      integer, intent(in) :: i, j, k, len
      integer, intent(in), dimension(3,(i+1)*(j+1)) :: pow

      real *8, intent(out), dimension((i+1)*(j+1)) :: coeffs_updt

      real *8, dimension(len) :: nck_coeffs
      integer, dimension(3,len) :: nck_pow
      integer, dimension(2,len) :: nck_der
      integer, dimension(len) :: nck_idx
      integer, dimension(2,k+1) :: nck_der_init
      integer, dimension(k+1) :: nck_idx_init
      integer :: l, l_idx
      
      logical, dimension(k+1) :: idx_flag0
      
      idx_flag0 = (2*(/(l, l=k,0,-1)/).le.i) .and.
     1            (2*(/(l, l=0,k, 1)/).le.j) 
      nck_der_init(1,:) = 2*(/(l, l=k,0,-1)/)
      nck_der_init(2,:) = 2*(/(l, l=0,k, 1)/)
      nck_idx_init = (/(l, l=1,k+1, 1)/)
      nck_der(1,:) = pack(nck_der_init(1,:),idx_flag0)
      nck_der(2,:) = pack(nck_der_init(2,:),idx_flag0)
      nck_idx = pack(nck_idx_init,idx_flag0)

      call polylapder2d(i, j, k, len, nck_der, nck_idx,
     1          nck_coeffs,nck_pow)
      coeffs_updt = 0
      do l=1,len
        do l_idx=1,(i+1)*(j+1)
          if ((pow(1,l_idx).eq.nck_pow(1,l)).and.
     1        (pow(2,l_idx).eq.nck_pow(2,l))) then ! power match, add contribution
            coeffs_updt(l_idx) = coeffs_updt(l_idx) + nck_coeffs(l)
          endif
        enddo
      enddo

      end subroutine hijcoeffsk0          


c     Fortran length of nck_coeffs fun. Hai 01/01/23
c
      subroutine lenofnckcoeffs0(i, j, ell, len)
c     
*/
      implicit none
      integer, intent(in) :: i, j, ell
      integer, intent(out) :: len
      
      logical, dimension(ell+1) :: idx_flag0
      integer :: k 
      
      idx_flag0 = (2*(/(k, k=ell,0, -1)/).le.i) .and.
     1            (2*(/(k, k=0,ell,  1)/).le.j) 
 
      len = count(idx_flag0)
      
      end subroutine lenofnckcoeffs0       

c     Fortran fun for the update part of Hij_coeffs from \Delta^k
c     Hai 01/01/23
c
      subroutine hijcoeffsk(i, j, k, pow, len, idx_flag, coeffs_updt)
c     
*/
      implicit none
      integer, intent(in) :: i, j, k, len
      integer, intent(in), dimension(2,(i+1)*(j+1)) :: pow
      integer, intent(in), dimension(k+1) :: idx_flag

      real *8, intent(out), dimension((i+1)*(j+1)) :: coeffs_updt

      real *8, dimension(len) :: nck_coeffs
      integer, dimension(3,len) :: nck_pow
      integer, dimension(2,len) :: nck_der
      integer, dimension(len) :: nck_idx
      integer, dimension(2,k+1) :: nck_der_init
      integer, dimension(k+1) :: nck_idx_init
      integer :: l, l_idx
      
      logical, dimension(k+1) :: idx_flag0
      
      idx_flag0 = idx_flag.eq.1 
      nck_der_init(1,:) = 2*(/(l, l=k,0, -1)/)
      nck_der_init(2,:) = 2*(/(l, l=0,k,  1)/)
      nck_idx_init = (/(l, l=1,k+1, 1)/)
      nck_der(1,:) = pack(nck_der_init(1,:),idx_flag0)
      nck_der(2,:) = pack(nck_der_init(2,:),idx_flag0)
      nck_idx = pack(nck_idx_init,idx_flag0)

      call polylapder2d(i, j, k, len, nck_der, nck_idx,
     1          nck_coeffs,nck_pow)
      coeffs_updt = 0
      do l=1,len
        do l_idx=1,(i+1)*(j+1)
          if ((pow(1,l_idx).eq.nck_pow(1,l)).and.
     1        (pow(2,l_idx).eq.nck_pow(2,l))) then ! power match, add contribution
            coeffs_updt(l_idx) = coeffs_updt(l_idx) + nck_coeffs(l)
          endif
        enddo
      enddo
      ! coeffs_updt(1) = nck_coeffs(1)
      ! coeffs_updt(2) = nck_coeffs(2)
      

      end subroutine hijcoeffsk          

c     Fortran length of nck_coeffs fun. Hai 01/01/23
c
      subroutine lenofnckcoeffs(i, j, ell, len, idx_flag)
c     
*/
      implicit none
      integer, intent(in) :: i, j, ell
      integer, intent(out) :: len
      integer, intent(out), dimension(ell+1) :: idx_flag

      logical, dimension(ell+1) :: idx_flag0
      integer :: k 
      
      idx_flag0 = (2*(/(k, k=ell,0, -1)/).le.i) .and.
     1            (2*(/(k, k=0,ell,  1)/).le.j) 
 
      len = count(idx_flag0)
      idx_flag = merge(1,0,idx_flag0)

      end subroutine lenofnckcoeffs        

c     Fortran polylapder2d fun. Hai 12/31/22. 
c     
      subroutine polylapder2d(i, j, ell, len, nck_der, nck_idx,
     1          nck_coeffs,nck_pow)
c     
*/
      implicit none
      integer, intent(in) :: i, j, ell, len
      integer, intent(in), dimension(2,len) :: nck_der
      integer, intent(in), dimension(len) :: nck_idx
      
      real *8, intent(out), dimension(len) :: nck_coeffs
      integer, intent(out), dimension(3,len) :: nck_pow

      integer :: k
      real *8 :: rec_coeffs

      do k=1,len
        nck_coeffs(k) = gamma(ell+1.0d0) / gamma(dble(nck_idx(k))) /
     1                  gamma(ell-nck_idx(k)+2.0d0)
        ! it is probably not that efficient to use gamma...
      enddo

      rec_coeffs = (-1)**ell*1.0d0/gamma(2.0d0*ell+2)
      do k=1,len
        nck_coeffs(k) = gamma(i+1.0d0)/gamma(i-nck_der(1,k)+1.0d0)
     1      *gamma(j+1.0d0)/gamma(j-nck_der(2,k)+1.0d0)*nck_coeffs(k)
      enddo
      nck_coeffs = rec_coeffs*nck_coeffs

      nck_pow(1,:) = i-nck_der(1,:)
      nck_pow(2,:) = j-nck_der(2,:)
      nck_pow(3,:) = (2*ell+1)
      

      end subroutine polylapder2d
      