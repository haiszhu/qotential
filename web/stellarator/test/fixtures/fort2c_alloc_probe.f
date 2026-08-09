      subroutine alloc_probe(n,ier)
      integer n,ier
      double precision, allocatable :: a(:)
      allocate(a(n),stat=ier)
      end
