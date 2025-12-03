      ! fortran call c... 
      ! to vectorize kernel, compute fsqrt?
      ! 07/15/24 Hai
      subroutine clap3ddlpmat(n,r0,m,r,rn,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: A(*)
      call clap3ddlpmat_c(n,r0,m,r,rn,w,A)  
      end subroutine clap3ddlpmat

      subroutine csimd128lap3ddlpmat(n,r0,m,r,rn,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: A(*)
      call csimd128lap3ddlpmat_c(n,r0,m,r,rn,w,A)  
      end subroutine csimd128lap3ddlpmat

      subroutine csimd256lap3ddlpmat(n,r0,m,r,rn,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: A(*)
      call csimd256lap3ddlpmat_c(n,r0,m,r,rn,w,A)   
      end subroutine csimd256lap3ddlpmat

      subroutine csimd512lap3ddlpmat(n,r0,m,r,rn,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: A(*)
      call csimd512lap3ddlpmat_c(n,r0,m,r,rn,w,A)
      end subroutine csimd512lap3ddlpmat

      subroutine clap3dslpmat(n,r0,m,r,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), w(*)
      real    *8, intent(inout) :: A(*)
      call clap3dslpmat_c(n,r0,m,r,w,A)  
      end subroutine clap3dslpmat

      subroutine csimd256lap3dslpmat(n,r0,m,r,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), w(*)
      real    *8, intent(inout) :: A(*)
      call csimd256lap3dslpmat_c(n,r0,m,r,w,A)
      end subroutine csimd256lap3dslpmat

      subroutine csimd512lap3dslpmat(n,r0,m,r,w,A)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), w(*)
      real    *8, intent(inout) :: A(*)
      call csimd512lap3dslpmat_c(n,r0,m,r,w,A)
      end subroutine csimd512lap3dslpmat

      subroutine csimd256lap3dsdlpmat(n,r0,m,r,rn,w,As,Ad)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: As(*),Ad(*)
      call csimd256lap3dsdlpmat_c(n,r0,m,r,rn,w,As,Ad)
      end subroutine csimd256lap3dsdlpmat

      subroutine csimd512lap3dsdlpmat(n,r0,m,r,rn,w,As,Ad)
      implicit none
      integer *8, intent(in) :: n, m
      real    *8, intent(in) :: r0(*), r(*), rn(*), w(*)
      real    *8, intent(inout) :: As(*),Ad(*)
      call csimd512lap3dsdlpmat_c(n,r0,m,r,rn,w,As,Ad)
      end subroutine csimd512lap3dsdlpmat

      subroutine cline3omega0123vr(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
      implicit none
      integer *8, intent(in) :: nbd, h_dim, order
      real    *8, intent(in) :: r0(*), m_all(*)
      real    *8, intent(in) :: onm0(*), onm1(*), onm2(*), onm3(*)
      real    *8, intent(inout) :: omega0(*),omega1(*)
      real    *8, intent(inout) :: omega2(*),omega3(*)
      call cline3omega0123vr_c(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
			end subroutine cline3omega0123vr

      subroutine csimd256line3omega0123vr(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
      implicit none
      integer *8, intent(in) :: nbd, h_dim, order
      real    *8, intent(in) :: r0(*), m_all(*)
      real    *8, intent(in) :: onm0(*), onm1(*), onm2(*), onm3(*)
      real    *8, intent(inout) :: omega0(*),omega1(*)
      real    *8, intent(inout) :: omega2(*),omega3(*)
      call csimd256line3omega0123vr_c(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
      end subroutine csimd256line3omega0123vr      

      subroutine csimd512line3omega0123vr(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
      implicit none
      integer *8, intent(in) :: nbd, h_dim, order
      real    *8, intent(in) :: r0(*), m_all(*)
      real    *8, intent(in) :: onm0(*), onm1(*), onm2(*), onm3(*)
      real    *8, intent(inout) :: omega0(*),omega1(*)
      real    *8, intent(inout) :: omega2(*),omega3(*)
      call csimd512line3omega0123vr_c(r0, m_all, nbd, order,    
     1                         onm0, onm1, onm2, onm3, h_dim, 
     2                         omega0, omega1, omega2, omega3)
      end subroutine csimd512line3omega0123vr

      subroutine cmomentsad_vr(r0, n, r, order, flag, Nk, Mk)
      implicit none
      integer *8, intent(in) :: n, order, flag
      real    *8, intent(in) :: r0(*), r(*)
      real    *8, intent(inout) :: Nk(*), Mk(*)
      call cmomentsad_vr_c(r0, n, r, order, flag, Nk, Mk)
      end subroutine cmomentsad_vr

      subroutine csimd512momentsad_vr(r0, n, r, order, flag, Nk, Mk)
      implicit none
      integer *8, intent(in) :: n, order, flag
      real    *8, intent(in) :: r0(*), r(*)
      real    *8, intent(inout) :: Nk(*), Mk(*)
      call csimd512momentsad_vr_c(r0, n, r, order, flag, Nk, Mk)
      end subroutine csimd512momentsad_vr

      subroutine csimd256momentsad_vrv2(r0, n, r, order, flag, Nk, Mk)
      implicit none
      integer *8, intent(in) :: n, order, flag
      real    *8, intent(in) :: r0(*), r(*)
      real    *8, intent(inout) :: Nk(*), Mk(*)
      call csimd256momentsad_vrv2_c(r0, n, r, order, flag, Nk, Mk)
      end subroutine csimd256momentsad_vrv2

      subroutine csimd512momentsad_vrv2(r0, n, r, order, flag, Nk, Mk)
      implicit none
      integer *8, intent(in) :: n, order, flag
      real    *8, intent(in) :: r0(*), r(*)
      real    *8, intent(inout) :: Nk(*), Mk(*)
      ! real    *8, intent(inout) :: Nk(n,order+1), Mk(n,order+1)
      ! integer *8 :: n8
      ! real    *8, allocatable :: Nk8(:,:), Mk8(:,:)
      call csimd512momentsad_vrv2_c(r0, n, r, order, flag, Nk, Mk)
      ! n8 = ((n+7)/8)*8
      ! allocate(Nk8(n8,order+1),Mk8(n8,order+1))
      ! call csimd512momentsad_vrv2_c(r0, n, r, order, flag, Nk8, Mk8)
      ! Nk = Nk8(1:n,:)
      ! Mk = Mk8(1:n,:)
      ! deallocate(Nk8,Mk8)
      end subroutine csimd512momentsad_vrv2

      ! subroutine lap3ddlpmat(n,r0,m,r,A)
      ! implicit none
      ! integer *8, intent(in) :: n, m
      ! real    *8, intent(in) :: r0(3,n), r(3,m)
      ! real    *8, intent(inout) :: A(n,m)
      ! A(1,1) = 0.0d0
      ! end subroutine lap3ddlpmat

      ! this is slower...
      subroutine frsqrt(x, rsqrt_x)
      use iso_fortran_env, only: int64
      implicit none
      real    *8, intent(in) :: x
      real    *8, intent(inout) :: rsqrt_x
      real    *8 :: x2, y
      integer *8, parameter :: ninter = 3
      integer *8 :: i
      real    *8, parameter :: half = 0.5d0
      real    *8, parameter :: onep5 = 1.5d0
      integer(int64), parameter :: magic = 6910469410427058089_int64
      x2 = half*x
      i  = transfer(x, i)
      i  = magic - shiftr(i, 1)
      y = transfer(i, rsqrt_x)
      do i = 1, ninter
        y = y * (onep5 - x2*y*y)  
      enddo
      rsqrt_x = y
      end subroutine frsqrt
