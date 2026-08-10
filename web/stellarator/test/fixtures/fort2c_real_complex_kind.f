      program main
      implicit none
      complex*16 z
      double precision x
      z = dcmplx(1.234567890123456d0,-2.345678901234567d0)
      x = real(z)
      print '(es24.16)', x
      end
