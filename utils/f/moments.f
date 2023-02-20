c     Fortran omega_all fun. Hai 01/06/23. compute vector part instead of scalar part...
c     other than that, it is exactly the same as omega_all
      subroutine omegaallvec(m, r0, dim1, m_all, n, morder, onm0, onm1,  
     1           onm2, onm3, h_dim, ijIdx, omegax, omegay, omegaz)
c     
*/    
      implicit none
      integer, intent(in) :: m, dim1, n, h_dim, morder
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(dim1,m) :: m_all
      real *8, intent(in), dimension(n*h_dim,4) :: onm0, onm1 
      real *8, intent(in), dimension(n*h_dim,4) :: onm2, onm3
      integer, intent(in), dimension(2,h_dim) :: ijIdx
      
      integer, dimension(h_dim) :: ijIdxsum
      integer :: j, k
      real *8, dimension(3) :: r0_j
      real *8, dimension(n,morder) :: mk_j
      real *8, dimension(n,h_dim) :: mkp1, mkp0
      real *8, dimension(n,h_dim) :: otmp
      real *8, dimension(h_dim) ::omega0_j,omega1_j,omega2_j,omega3_j
      real *8 :: pi, tmp
      real *8, intent(out), dimension(m,4*h_dim) :: omegax,omegay,omegaz
      integer, dimension(h_dim) :: tmp_vec
      real *8, dimension(n,h_dim,4) :: onm0_rs, onm1_rs
      real *8, dimension(n,h_dim,4) :: onm2_rs,onm3_rs
      
      pi = 4.0d0*atan(1.0d0)
      tmp = -1.0d0/(4.0d0*pi)
      tmp_vec = (/(k, k=1,h_dim, 1)/)

      ijIdxsum = sum(ijIdx,dim=1)
      onm0_rs = reshape(onm0, (/n, h_dim, 4/))
      onm1_rs = reshape(onm1, (/n, h_dim, 4/))
      onm2_rs = reshape(onm2, (/n, h_dim, 4/))
      onm3_rs = reshape(onm3, (/n, h_dim, 4/))

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(j,r0_j,mk_j,mkp1,mkp0,otmp)
C$OMP$PRIVATE(omega0_j,omega1_j,omega2_j,omega3_j)
      do j=1,m
        r0_j = r0(:,j)
        mk_j = reshape(M_all(:,j), (/n, morder/))
        
        mkp1 = mk_j(:,ijIdxsum+3)
        mkp0 = mk_j(:,ijIdxsum+2)

        otmp = onm0_rs(:,:,2)*r0_j(1) + onm0_rs(:,:,3)*r0_j(2) 
     1         + onm0_rs(:,:,4)*r0_j(3)
        omega0_j= tmp*sum(mkp1*onm0_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm1_rs(:,:,2)*r0_j(1) + onm1_rs(:,:,3)*r0_j(2)
     1         + onm1_rs(:,:,4)*r0_j(3)
        omega1_j = tmp*sum(mkp1*onm1_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm2_rs(:,:,2)*r0_j(1) + onm2_rs(:,:,3)*r0_j(2) 
     1         + onm2_rs(:,:,4)*r0_j(3)
        omega2_j = tmp*sum(mkp1*onm2_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm3_rs(:,:,2)*r0_j(1) + onm3_rs(:,:,3)*r0_j(2)
     1         + onm3_rs(:,:,4)*r0_j(3)
        omega3_j = tmp*sum(mkp1*onm3_rs(:,:,1)+mkp0*otmp,dim=1)

        omegax(j,tmp_vec) = omega1_j
        omegax(j,tmp_vec+h_dim) = omega0_j
        omegax(j,tmp_vec+2*h_dim) = -omega3_j
        omegax(j,tmp_vec+3*h_dim) = omega2_j

        omegay(j,tmp_vec) = omega2_j
        omegay(j,tmp_vec+h_dim) = omega3_j
        omegay(j,tmp_vec+2*h_dim) = omega0_j
        omegay(j,tmp_vec+3*h_dim) = -omega1_j

        omegaz(j,tmp_vec) = omega3_j
        omegaz(j,tmp_vec+h_dim) = -omega2_j
        omegaz(j,tmp_vec+2*h_dim) = omega1_j
        omegaz(j,tmp_vec+3*h_dim) = omega0_j

      enddo
C$OMP END PARALLEL DO   

      end subroutine omegaallvec

c     Fortran omega_all fun. Hai 12/29/22. should be put into omega later?
      subroutine omegaall(m, r0, dim1, m_all, n, morder, onm0, onm1,  
     1           onm2, onm3, h_dim, ijIdx, omega)
c     
*/    
      implicit none
      integer, intent(in) :: m, dim1, n, h_dim, morder
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(dim1,m) :: m_all
      real *8, intent(in), dimension(n*h_dim,4) :: onm0, onm1 
      real *8, intent(in), dimension(n*h_dim,4) :: onm2, onm3
      integer, intent(in), dimension(2,h_dim) :: ijIdx
      
      integer, dimension(h_dim) :: ijIdxsum
      integer :: j, k
      real *8, dimension(3) :: r0_j
      real *8, dimension(n,morder) :: mk_j
      real *8, dimension(n,h_dim) :: mkp1, mkp0
      real *8, dimension(n,h_dim) :: otmp
      real *8, dimension(h_dim) ::omega0_j,omega1_j,omega2_j,omega3_j
      real *8 :: pi, tmp
      real *8, intent(out), dimension(m,4*h_dim) :: omega
      integer, dimension(h_dim) :: tmp_vec
      real *8, dimension(n,h_dim,4) :: onm0_rs, onm1_rs
      real *8, dimension(n,h_dim,4) :: onm2_rs,onm3_rs
      
      pi = 4.0d0*atan(1.0d0)
      tmp = -1.0d0/(4.0d0*pi)
      tmp_vec = (/(k, k=1,h_dim, 1)/)

      ijIdxsum = sum(ijIdx,dim=1)
      onm0_rs = reshape(onm0, (/n, h_dim, 4/))
      onm1_rs = reshape(onm1, (/n, h_dim, 4/))
      onm2_rs = reshape(onm2, (/n, h_dim, 4/))
      onm3_rs = reshape(onm3, (/n, h_dim, 4/))

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(j,r0_j,mk_j,mkp1,mkp0,otmp)
C$OMP$PRIVATE(omega0_j,omega1_j,omega2_j,omega3_j)
      do j=1,m
        r0_j = r0(:,j)
        mk_j = reshape(M_all(:,j), (/n, morder/))
        
        mkp1 = mk_j(:,ijIdxsum+3)
        mkp0 = mk_j(:,ijIdxsum+2)

        otmp = onm0_rs(:,:,2)*r0_j(1) + onm0_rs(:,:,3)*r0_j(2) 
     1         + onm0_rs(:,:,4)*r0_j(3)
        omega0_j= tmp*sum(mkp1*onm0_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm1_rs(:,:,2)*r0_j(1) + onm1_rs(:,:,3)*r0_j(2)
     1         + onm1_rs(:,:,4)*r0_j(3)
        omega1_j = tmp*sum(mkp1*onm1_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm2_rs(:,:,2)*r0_j(1) + onm2_rs(:,:,3)*r0_j(2) 
     1         + onm2_rs(:,:,4)*r0_j(3)
        omega2_j = tmp*sum(mkp1*onm2_rs(:,:,1)+mkp0*otmp,dim=1)

        otmp = onm3_rs(:,:,2)*r0_j(1) + onm3_rs(:,:,3)*r0_j(2)
     1         + onm3_rs(:,:,4)*r0_j(3)
        omega3_j = tmp*sum(mkp1*onm3_rs(:,:,1)+mkp0*otmp,dim=1)

        omega(j,tmp_vec) = omega0_j
        omega(j,tmp_vec+h_dim) = -omega1_j
        omega(j,tmp_vec+2*h_dim) = -omega2_j
        omega(j,tmp_vec+3*h_dim) = -omega3_j

      enddo
C$OMP END PARALLEL DO   

      end subroutine omegaall

c     Fortran moments_all fun, non adaptive version. Hai 02/20/23
      subroutine momentsallplain(m, r0, n, r, order, np, p, dim1,
     1         m_all)
c     
*/    
      implicit none
      integer, intent(in) :: m, n, order, np, p, dim1
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(3,n) :: r
      real *8, intent(out), dimension(dim1,m) :: m_all

      integer :: j, ell, k
      integer, dimension(p) :: idx_ell, tmp_vec
      real *8, dimension(3) :: r0_j
      real *8, dimension(3,p) :: r_ell
      real *8, dimension(p,3) :: r_ell_T
      real *8, dimension(p,2*order+2) :: nk_j, mk_j

      tmp_vec = (/(k, k=1,p, 1)/) 

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(j,ell,k,r0_j,idx_ell,r_ell,r_ell_T,nk_j,mk_j)
      do j=1,m
        r0_j = r0(:,j)
        do ell=1,np
          idx_ell = (ell-1)*p + tmp_vec
          r_ell = r(:,idx_ell)
          r_ell_T = transpose(r_ell)
          call moments(r0_j, p, r_ell_T, 2*order+1, 0, nk_j, mk_j)
          do k=1,2*order+2
            m_all((k-1)*n+idx_ell,j) = mk_j(:,k)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO   

      end subroutine momentsallplain  

c     Fortran moments_all fun. Hai 12/29/22
c     add later

c     Fortran moments_up fun. Hai 12/29/22
c     add later

c     part of Fortran moments_up fun. Hai 12/28/22
c     add later

c     Fortran line3_adaptive_r_near_r0 fun. Hai 12/28/22
c     add later

c     Fortran length of pan_t_end
c     add later

c     Fortran add_extra_panel fun. Hai 12/28/22
c     add later

c     Fortran bclag_interp_matrix fun. Hai 12/27/22
c     add later


c     Fortran moments fun. Hai 12/27/22

      subroutine moments(r0, n, r, order, flag, Nk, Mk)
c     
*/
      implicit none
      integer, intent(in) :: n, flag, order
      real *8, intent(in), dimension(3) :: r0
      real *8, intent(in), dimension(n,3) :: r
      real *8, intent(out), dimension(n,order+1) :: Nk, Mk

      integer :: k
      real *8 :: r0norm, r0norm_inv
      real *8, dimension(n) :: r0dotr, rnorm, rnorm_inv, rnorm2_inv
      real *8, dimension(n) :: r0dotr_over_rnorm2, r0norm2_over_rnorm2
      real *8, dimension(n,3) :: r0mr_vec
      real *8, dimension(n) :: r0mr, r0mr_inv, r0mr_over_rnorm2
      real *8, dimension(n) :: r0dotr0mr_over_r0mr, rdotr0mr_over_r0mr
      real *8, dimension(n) :: denominator1, denominator2, LMNcommon
      real *8, dimension(n) :: N0, N1, M0, M1
      real *8, dimension(n) :: LMcommon

      real *8, dimension(n) :: r0normplusrnorm
      
      r0dotr = r0(1)*r(:,1) + r0(2)*r(:,2) + r0(3)*r(:,3)
      rnorm = sqrt(sum(r**2, dim=2))
      r0norm = sqrt(r0(1)**2 + r0(2)**2 + r0(3)**2)
      rnorm_inv = 1/rnorm
      rnorm2_inv = rnorm_inv**2
      r0norm_inv = 1/r0norm
      r0dotr_over_rnorm2 = r0dotr*rnorm2_inv
      r0norm2_over_rnorm2 = r0norm**2*rnorm2_inv
      r0normplusrnorm = r0norm+rnorm

      r0mr_vec(:,1) = r0(1) - r(:,1)
      r0mr_vec(:,2) = r0(2) - r(:,2)
      r0mr_vec(:,3) = r0(3) - r(:,3)
      r0mr = sqrt(sum(r0mr_vec**2, dim=2))
      r0mr_over_rnorm2 = r0mr*rnorm2_inv;
      r0mr_inv = 1/r0mr;

c     shared by N, M, & L   
      r0dotr0mr_over_r0mr = (r0(1)*r0mr_vec(:,1)+r0(2)*r0mr_vec(:,2)+ 
     1                       r0(3)*r0mr_vec(:,3))*r0mr_inv
      denominator1 = r0norm + r0dotr0mr_over_r0mr
      rdotr0mr_over_r0mr = (r(:,1)*r0mr_vec(:,1)+r(:,2)*r0mr_vec(:,2)+ 
     1                       r(:,3)*r0mr_vec(:,3))*r0mr_inv
      denominator2 = rnorm + rdotr0mr_over_r0mr
      LMNcommon = r0normplusrnorm/(denominator1+denominator2)

cc    N
c     \int t^0/|r0-t*r| dt
      N0 = log((r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv)*rnorm_inv
      ! N0 = (-log(r0mr)+log((r0norm+rnorm+r0mr)*LMNcommon))*rnorm_inv
c     \int t^1/|r0-t*r| dt
      N1 = N0*r0dotr_over_rnorm2 + (r0mr-r0norm)*rnorm2_inv
c     initialize \int t^k/|r0-t*r| dt
      Nk(:,1) = N0
      Nk(:,2) = N1
c     recursion 
      do k=2,order
        Nk(:,k+1) = dble(2*k-1)/k*r0dotr_over_rnorm2*Nk(:,k) - 
     1              dble(k-1)/k*r0norm2_over_rnorm2*Nk(:,k-1) + 
     2              dble(1)/k*r0mr_over_rnorm2           
      enddo

cc    M
      LMcommon = 1/(r0norm*rnorm+r0dotr)
c     \int t^0/|r0-t*r|^3 dt
      M0 = LMNcommon*LMcommon*(
     1((r0normplusrnorm)*r0mr_inv+rnorm*r0norm_inv)*r0mr_inv-r0norm_inv)
c     \int t^1/|r0-t*r|^3 dt
      M1 = r0norm*M0/(r0norm+r0mr)
c     initialize \int t^k/|r0-t*r|^3 dt
      Mk(:,1) = M0 
      Mk(:,2) = M1
c     recursion
      do k=2,order
        Mk(:,k+1) = (r0dotr*Mk(:,k)+(k-1)*Nk(:,k-1)-r0mr_inv)*rnorm2_inv
      enddo
      
      end subroutine moments

ccccccccccccccccc below are for Laplace DLPn, should be combined to code above...

c     Fortran moments_all fun. Hai 01/14/23
c     add later


c     Fortran moments_up fun. Hai 12/29/22
c     add later


c     Fortran moments fun with Lk for DLPn. Hai 01/14/23

      subroutine moments2(r0, n, r, order, flag, Nk, Mk, Lk)
c     
*/
      implicit none
      integer, intent(in) :: n, flag, order
      real *8, intent(in), dimension(3) :: r0
      real *8, intent(in), dimension(n,3) :: r
      real *8, intent(out), dimension(n,order+1) :: Nk, Mk, Lk

      integer :: k
      real *8 :: r0norm, r0norm_inv
      real *8, dimension(n) :: r0dotr, rnorm, rnorm_inv, rnorm2_inv
      real *8, dimension(n) :: r0dotr_over_rnorm2, r0norm2_over_rnorm2
      real *8, dimension(n,3) :: r0mr_vec
      real *8, dimension(n) :: r0mr, r0mr_inv, r0mr_over_rnorm2
      real *8, dimension(n) :: r0dotr0mr_over_r0mr, rdotr0mr_over_r0mr
      real *8, dimension(n) :: denominator1, denominator2, LMNcommon
      real *8, dimension(n) :: N0, N1, M0, M1, L0, L1
      real *8, dimension(n) :: LMcommon, Lcommon1, Lcommon2

      real *8, dimension(n) :: r0normplusrnorm
      
      r0dotr = r0(1)*r(:,1) + r0(2)*r(:,2) + r0(3)*r(:,3)
      rnorm = sqrt(sum(r**2, dim=2))
      r0norm = sqrt(r0(1)**2 + r0(2)**2 + r0(3)**2)
      rnorm_inv = 1/rnorm
      rnorm2_inv = rnorm_inv**2
      r0norm_inv = 1/r0norm
      r0dotr_over_rnorm2 = r0dotr*rnorm2_inv
      r0norm2_over_rnorm2 = r0norm**2*rnorm2_inv
      r0normplusrnorm = r0norm+rnorm

      r0mr_vec(:,1) = r0(1) - r(:,1)
      r0mr_vec(:,2) = r0(2) - r(:,2)
      r0mr_vec(:,3) = r0(3) - r(:,3)
      r0mr = sqrt(sum(r0mr_vec**2, dim=2))
      r0mr_over_rnorm2 = r0mr*rnorm2_inv;
      r0mr_inv = 1/r0mr;

c     shared by N, M, & L   
      r0dotr0mr_over_r0mr = (r0(1)*r0mr_vec(:,1)+r0(2)*r0mr_vec(:,2)+ 
     1                       r0(3)*r0mr_vec(:,3))*r0mr_inv
      denominator1 = r0norm + r0dotr0mr_over_r0mr
      rdotr0mr_over_r0mr = (r(:,1)*r0mr_vec(:,1)+r(:,2)*r0mr_vec(:,2)+ 
     1                       r(:,3)*r0mr_vec(:,3))*r0mr_inv
      denominator2 = rnorm + rdotr0mr_over_r0mr
      LMNcommon = r0normplusrnorm/(denominator1+denominator2)

cc    N
c     \int t^0/|r0-t*r| dt
      N0 = log((r0normplusrnorm+r0mr)*LMNcommon*r0mr_inv)*rnorm_inv
      ! N0 = (-log(r0mr)+log((r0norm+rnorm+r0mr)*LMNcommon))*rnorm_inv
c     \int t^1/|r0-t*r| dt
      N1 = N0*r0dotr_over_rnorm2 + (r0mr-r0norm)*rnorm2_inv
c     initialize \int t^k/|r0-t*r| dt
      Nk(:,1) = N0
      Nk(:,2) = N1
c     recursion 
      do k=2,order
        Nk(:,k+1) = dble(2*k-1)/k*r0dotr_over_rnorm2*Nk(:,k) - 
     1              dble(k-1)/k*r0norm2_over_rnorm2*Nk(:,k-1) + 
     2              dble(1)/k*r0mr_over_rnorm2           
      enddo

cc    M
      LMcommon = 1/(r0norm*rnorm+r0dotr)
c     \int t^0/|r0-t*r|^3 dt
      M0 = LMNcommon*LMcommon*(
     1((r0normplusrnorm)*r0mr_inv+rnorm*r0norm_inv)*r0mr_inv-r0norm_inv)
c     \int t^1/|r0-t*r|^3 dt
      M1 = r0norm*M0/(r0norm+r0mr)
c     initialize \int t^k/|r0-t*r|^3 dt
      Mk(:,1) = M0 
      Mk(:,2) = M1
c     recursion
      do k=2,order
        Mk(:,k+1) = (r0dotr*Mk(:,k)+(k-1)*Nk(:,k-1)-r0mr_inv)*rnorm2_inv
      enddo

cc    L
      Lcommon1 = r0norm**2*LMNcommon**2*(rnorm+rdotr0mr_over_r0mr/3
     1            -r0dotr*r0norm_inv/3)
      Lcommon2 = r0norm*LMNcommon*(-rnorm+rdotr0mr_over_r0mr/3
     1            +2*r0dotr*r0norm_inv/3)
c     \int t^0/|r0-t*r|^5 dt
      L0=((((r0norm+rnorm)*Lcommon1*((r0norm+rnorm)*r0mr_inv+2)*r0mr_inv 
     1 +(Lcommon1+(r0norm+rnorm)*Lcommon2))*r0mr_inv+Lcommon2)*r0mr_inv 
     2 -(2*rdotr0mr_over_r0mr/3+r0dotr*r0norm_inv/3))*LMcommon**2
     3  *r0norm_inv**2
c     \int t^1/|r0-t*r|^5 dt
      L1 = (L0*r0dotr - (r0mr_inv**3-r0norm_inv**3)/3)*rnorm2_inv
c     initialize \int t^k/|r0-t*r|^5 dt
      Lk(:,1) = L0 
      Lk(:,2) = L1
c     recursion
      do k=2,order
        Lk(:,k+1) = (2*r0dotr*Lk(:,k) - r0norm**2*Lk(:,k-1) 
     1               + Mk(:,k-1))*rnorm2_inv
      enddo

      end subroutine moments2

c     Fortran omega for DLPn fun. Hai 01/14/23.
c     add later


c     Fortran omega for all 9 components of DLPn fun. Hai 01/18/23.
c     add later
