c     Lap3dSLPmat... Hai 01/17/23
c
      subroutine lap3dslpmat(m,r0,n,r,w,A)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(3,n) :: r
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n) :: A
      
      integer :: j,k
      real *8 :: pi
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(m) :: r0x, r0y, r0z

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = w/(4.0d0*pi)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
      do k=1,n
        do j=1,m
          A(j,k) = wtmp(k)/sqrt
     1    ((r0x(j)-rx(k))**2+(r0y(j)-ry(k))**2+(r0z(j)-rz(k))**2)
        enddo
      enddo 
C$OMP END PARALLEL DO   
      end subroutine lap3dslpmat

      subroutine lap3ddlpmat(m,r0,n,r,rn,w,A)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(3,n) :: r, rn
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n) :: A

      integer :: j,k
      real *8 :: pi, rr, dx, dy, dz
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(n) :: rnx, rny, rnz
      real *8, dimension(m) :: r0x, r0y, r0z

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)
      rnx = rn(1,:) 
      rny = rn(2,:) 
      rnz = rn(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = w/(4.0d0*pi)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
C$OMP$PRIVATE(dx,dy,dz,rr)
      do k=1,n
        do j=1,m
          dx = r0x(j)-rx(k)
          dy = r0y(j)-ry(k)
          dz = r0z(j)-rz(k)
          rr = (dx**2+dy**2+dz**2)
          A(j,k) = wtmp(k)*(dx*rnx(k)+dy*rny(k)+dz*rnz(k))/(sqrt(rr)*rr) 
        enddo
      enddo 
C$OMP END PARALLEL DO     

      end subroutine lap3ddlpmat


      subroutine sto3ddlpmat(m,r0,n,r,rn,w,A11,A12,A13,A22,A23,A33)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(3,n) :: r, rn
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n)::A11, A12, A13, A22, A23, A33

      integer :: j,k
      real *8 :: pi, ir, ir2, dx, dy, dz, dnxir2
      real *8 :: dxdx, dydy, dzdz
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(n) :: rnx, rny, rnz
      real *8, dimension(m) :: r0x, r0y, r0z

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)
      rnx = rn(1,:) 
      rny = rn(2,:) 
      rnz = rn(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = 3.0d0*w/(4.0d0*pi)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
C$OMP$PRIVATE(dx,dy,dz,dxdx,dydy,dzdz,ir,ir2,dnxir2)
      do k=1,n
        do j=1,m
          dx = r0x(j)-rx(k)
          dy = r0y(j)-ry(k)
          dz = r0z(j)-rz(k)
          dxdx = dx*dx
          dydy = dy*dy
          dzdz = dz*dz 
          ir = 1.0d0/sqrt(dxdx+dydy+dzdz)
          ir2 = ir**2
          dnxir2 = (dx*rnx(k)+dy*rny(k)+dz*rnz(k))*ir2*ir*ir2*wtmp(k)

          A11(j,k) = dxdx*dnxir2
          A12(j,k) = dx*dy*dnxir2
          A13(j,k) = dx*dz*dnxir2
          A22(j,k) = dydy*dnxir2
          A23(j,k) = dy*dz*dnxir2
          A33(j,k) = dzdz*dnxir2
        enddo
      enddo   
C$OMP END PARALLEL DO  

      end subroutine sto3ddlpmat

      subroutine sto3ddlpnmat(m,r0,r0n,n,r,rn,w,A11,A12,A13,A22,A23,A33,
     1                        A21,A31,A32)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0, r0n
      real *8, intent(in), dimension(3,n) :: r, rn
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n)::A11, A12, A13, A22, A23, A33
      real *8, intent(out), dimension(m,n)::A21, A31, A32

      integer :: j,k
      real *8 :: pi, ir, ir2, ir3, ir5, dx, dy, dz, dnxir2, dnx, Anjk
      real *8 :: dnxir5, dnyir5, dnxdnyir5, nxnyir5
      real *8 :: nyir31, nyir32, nyir33
      real *8 :: ddnxir51, ddnxir52, ddnxir53
      real *8 :: ddnyir51, ddnyir52, ddnyir53
      real *8 :: dxdx, dydy, dzdz
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(n) :: rnx, rny, rnz
      real *8, dimension(m) :: r0x, r0y, r0z
      real *8, dimension(m) :: r0nx, r0ny, r0nz

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      r0nx = r0n(1,:) 
      r0ny = r0n(2,:) 
      r0nz = r0n(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)
      rnx = rn(1,:) 
      rny = rn(2,:) 
      rnz = rn(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = 3.0d0*w/(4.0d0*pi)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
C$OMP$PRIVATE(Anjk,dx,dy,dz,dxdx,dydy,dzdz,ir,ir2,ir3,ir5,dnxir2)
C$OMP$PRIVATE(dnx,dnxir5,dnyir5,dnxdnyir5,nxnyir5)
C$OMP$PRIVATE(nyir31,nyir32,nyir33)
C$OMP$PRIVATE(ddnxir51,ddnxir52,ddnxir53)
C$OMP$PRIVATE(ddnyir51,ddnyir52,ddnyir53)
      do k=1,n
        do j=1,m
          dx = r0x(j)-rx(k)
          dy = r0y(j)-ry(k)
          dz = r0z(j)-rz(k)
          dxdx = dx*dx
          dydy = dy*dy
          dzdz = dz*dz 
          ir = 1.0d0/sqrt(dxdx+dydy+dzdz)
          ir2 = ir**2
          dnxir2 = (dx*rnx(k)+dy*rny(k)+dz*rnz(k))*ir2*ir*ir2*wtmp(k)

          ir3 = ir2*ir*w(k)
          ir5 = ir2**2*ir*wtmp(k)
          dnx = (dx*r0nx(j)+dy*r0ny(j)+dz*r0nz(j))
          dnxir5 = dnx*ir5
          dnyir5 = (dx*rnx(k)+dy*rny(k)+dz*rnz(k))*ir5

          nxnyir5 = (r0nx(j)*rnx(k)+r0ny(j)*rny(k)+r0nz(j)*rnz(k))*ir5
          dnxdnyir5 = dnx*dnyir5

          Anjk = nxnyir5-10*dnxdnyir5*ir2

          nyir31 = rnx(k)*ir3
          nyir32 = rny(k)*ir3
          nyir33 = rnz(k)*ir3

          ddnxir51 = dx*dnxir5*(4.0d0*pi)/(3.0d0)
          ddnxir52 = dy*dnxir5*(4.0d0*pi)/(3.0d0)
          ddnxir53 = dz*dnxir5*(4.0d0*pi)/(3.0d0)

          ddnyir51 = dx*dnyir5*(4.0d0*pi)/(3.0d0)
          ddnyir52 = dy*dnyir5*(4.0d0*pi)/(3.0d0)
          ddnyir53 = dz*dnyir5*(4.0d0*pi)/(3.0d0)

          A11(j,k) =  2*r0nx(j)*nyir31 + 3*rnx(k)*ddnxir51 
     1              + 3*ddnyir51*r0nx(j)
          A12(j,k) =  2*r0nx(j)*nyir32 + 3*rnx(k)*ddnxir52 
     1              + 3*ddnyir51*r0ny(j)
          A13(j,k) =  2*r0nx(j)*nyir33 + 3*rnx(k)*ddnxir53 
     1              + 3*ddnyir51*r0nz(j)
          A22(j,k) =  2*r0ny(j)*nyir32 + 3*rny(k)*ddnxir52 
     1              + 3*ddnyir52*r0ny(j)
          A23(j,k) =  2*r0ny(j)*nyir33 + 3*rny(k)*ddnxir53 
     1              + 3*ddnyir52*r0nz(j)
          A33(j,k) =  2*r0nz(j)*nyir33 + 3*rnz(k)*ddnxir53 
     1              + 3*ddnyir53*r0nz(j)


          A11(j,k) = 1.0d0/(4.0d0*pi)*A11(j,k) + Anjk*dxdx + dnxdnyir5
          A12(j,k) = 1.0d0/(4.0d0*pi)*A12(j,k) + Anjk*dx*dy
          A13(j,k) = 1.0d0/(4.0d0*pi)*A13(j,k) + Anjk*dx*dz
          A22(j,k) = 1.0d0/(4.0d0*pi)*A22(j,k) + Anjk*dydy + dnxdnyir5
          A23(j,k) = 1.0d0/(4.0d0*pi)*A23(j,k) + Anjk*dy*dz
          A33(j,k) = 1.0d0/(4.0d0*pi)*A33(j,k) + Anjk*dzdz + dnxdnyir5

          A21(j,k) =  2*r0ny(j)*nyir31 + 3*rny(k)*ddnxir51 
     1              + 3*ddnyir52*r0nx(j)
          A31(j,k) =  2*r0nz(j)*nyir31 + 3*rnz(k)*ddnxir51 
     1              + 3*ddnyir53*r0nx(j)
          A32(j,k) =  2*r0nz(j)*nyir32 + 3*rnz(k)*ddnxir52 
     1              + 3*ddnyir53*r0ny(j)
          
          A21(j,k) = 1.0d0/(4.0d0*pi)*A21(j,k) + Anjk*dx*dy
          A31(j,k) = 1.0d0/(4.0d0*pi)*A31(j,k) + Anjk*dx*dz
          A32(j,k) = 1.0d0/(4.0d0*pi)*A32(j,k) + Anjk*dy*dz

        enddo
      enddo   
C$OMP END PARALLEL DO  

      end subroutine sto3ddlpnmat

      subroutine sto3dslpmat(m,r0,n,r,w,A11,A12,A13,A22,A23,A33)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0
      real *8, intent(in), dimension(3,n) :: r
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n)::A11, A12, A13, A22, A23, A33

      integer :: j,k
      real *8 :: pi, ir, ir3, dx, dy, dz
      real *8 :: dxdx, dydy, dzdz
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(m) :: r0x, r0y, r0z

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = 1.0d0*w/(8.0d0*pi)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
C$OMP$PRIVATE(dx,dy,dz,dxdx,dydy,dzdz,ir,ir3)
      do k=1,n
        do j=1,m
          dx = r0x(j)-rx(k)
          dy = r0y(j)-ry(k)
          dz = r0z(j)-rz(k)
          dxdx = dx*dx
          dydy = dy*dy
          dzdz = dz*dz 
          ir = 1.0d0/sqrt(dxdx+dydy+dzdz)
          ir3 = ir**3*wtmp(k)
          ir = ir*wtmp(k)

          A11(j,k) = dxdx*ir3+ir
          A12(j,k) = dx*dy*ir3
          A13(j,k) = dx*dz*ir3
          A22(j,k) = dydy*ir3+ir
          A23(j,k) = dy*dz*ir3
          A33(j,k) = dzdz*ir3+ir
        enddo
      enddo   
C$OMP END PARALLEL DO  

      end subroutine sto3dslpmat


      subroutine sto3dslpnmat(m,r0,r0n,n,r,w,A11,A12,A13,A22,A23,A33)
c     
*/    
      implicit none
      integer, intent(in) :: m, n
      real *8, intent(in), dimension(3,m) :: r0, r0n
      real *8, intent(in), dimension(3,n) :: r
      real *8, intent(in), dimension(n) :: w
      real *8, intent(out), dimension(m,n)::A11, A12, A13, A22, A23, A33

      integer :: j,k
      real *8 :: pi, ir, ir2, dx, dy, dz, dnxir2
      real *8 :: dxdx, dydy, dzdz
      real *8, dimension(n) :: wtmp, rx, ry, rz
      real *8, dimension(m) :: r0nx, r0ny, r0nz
      real *8, dimension(m) :: r0x, r0y, r0z

      r0x = r0(1,:)
      r0y = r0(2,:)
      r0z = r0(3,:)
      rx = r(1,:)
      ry = r(2,:)
      rz = r(3,:)
      r0nx = r0n(1,:) 
      r0ny = r0n(2,:) 
      r0nz = r0n(3,:)

      pi = 4.0d0*atan(1.0d0)

      wtmp = -3.0d0*w/(4.0d0*pi)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(k,j)
C$OMP$PRIVATE(dx,dy,dz,dxdx,dydy,dzdz,ir,ir2,dnxir2)
      do k=1,n
        do j=1,m
          dx = r0x(j)-rx(k)
          dy = r0y(j)-ry(k)
          dz = r0z(j)-rz(k)
          dxdx = dx*dx
          dydy = dy*dy
          dzdz = dz*dz 
          ir = 1.0d0/sqrt(dxdx+dydy+dzdz)
          ir2 = ir**2
          dnxir2 = (dx*r0nx(j)+dy*r0ny(j)+dz*r0nz(j))*ir2*ir*ir2*wtmp(k)

          A11(j,k) = dxdx*dnxir2
          A12(j,k) = dx*dy*dnxir2
          A13(j,k) = dx*dz*dnxir2
          A22(j,k) = dydy*dnxir2
          A23(j,k) = dy*dz*dnxir2
          A33(j,k) = dzdz*dnxir2
        enddo
      enddo   
C$OMP END PARALLEL DO  

      end subroutine sto3dslpnmat