module module_math
  use module_globals
  use module_autodiff
  implicit none

  private
  public IntegralTest, load_integral_test, f, Df

  type :: IntegralTest
     real(8) :: st, en, le, cn, seikai
  end type IntegralTest
  
contains
  !========================================================================================
  subroutine load_integral_test(test)
    type(IntegralTest),intent(out) :: test

    ! int_[a->b](sin(x))dx
    ! test%st=0.d0
    ! test%en=30.d0
    ! test%seikai=-cos(test%en)+cos(test%st)

    !int_[0.01->2](sin(x)/x)dx
    ! test%st=0.01d0
    ! test%en=1.d0
    ! test%seikai=0.93608312592257190411368860853800106052127274d0

    ! !int_[0->30](sin(x)/x)dx
    ! test%st=1.d0
    ! test%en=30.d0
    ! test%seikai=0.620673469663168096042377995183618508711157d0

    ! int_[-1->1] 1/(1+25*x^2) dx
    test%st=-1.d0
    test%en=1.d0
    test%seikai=(atan(5.d0*test%en)-atan(5.d0*test%st))/5.d0

    ! int_[0,1] bessel_j0(x) dx
    ! test%st=0.d0
    ! test%en=1.d0
    ! test%seikai=0.919730410089760239314421194080619970661d0

    ! test%st=0.d0
    ! test%en=2.d0
    ! test%seikai=(test%en**30-test%st**30)/30.d0
    
    test%le=test%en-test%st
    test%cn=(test%st+test%en)*0.5d0

  end subroutine load_integral_test

  !========================================================================================
  real(8) function f(x)
    real(8),intent(in) :: x
    !f=sin(x)
    !f=sin(x)/x
    !f=sin(x)/x
    f=1.d0/(1.d0+25.d0*x*x)
    !f=bessel_j0(x)
    !f=x**29
  end function f

  !========================================================================================
  type(Dcomplex8) function Df(x)
    type(Dcomplex8),intent(in) :: x
    !Df=sin(x)
    !Df=sin(x)/x
    !Df=sin(x)/x
    Df=1.d0/(1.d0+25.d0*x*x)
    !Df=bessel_j0(x)
    !Df=x**29
  end function Df
  
end module module_math
