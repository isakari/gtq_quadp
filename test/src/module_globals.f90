module module_globals
  implicit none

  complex(8),parameter :: one=cmplx(1.d0,0.d0,kind(1.d0))
  complex(8),parameter :: zero=cmplx(0.d0,0.d0,kind(1.d0))
  complex(8),parameter :: ione=cmplx(0.d0,1.d0,kind(1.d0))
  real(8),parameter :: pi=acos(-1.d0)
  integer,parameter :: ntaylor=20
  integer,parameter :: maxs=ntaylor/2
  integer,parameter :: maxngt=10
  
end module module_globals
