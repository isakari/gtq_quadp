program main
  use module_globals
  use module_autodiff
  use module_math
  use module_gauss_turan_quadrature
  implicit none
  
  type(GaussTuranQuadrature) :: gt(maxngt)
  type(iCj) :: choose(0:ntaylor)
  type(IntegralTest) :: test

  integer :: ng, i, is, k
  real(8) :: og, coef

  type(Dcomplex8) :: ff, xx

  call load_integral_test(test)
  call init_gauss_turan(maxs,maxngt,gt)
  call init_autodiff(choose)

  do is=0,maxs
     do ng=1,maxngt
        og=0.d0
        do i=1,ng
           xx%f(0)=cmplx(test%le*0.5d0*gt(ng)%s(is)%gz(i)+test%cn,0.d0,kind(1.d0))
           xx%f(1)=one
           xx%f(2:ntaylor)=zero
           ff=Df(xx)
           coef=0.5d0*test%le
           do k=0,2*is
              og=og+real(ff%f(k))*gt(ng)%s(is)%we(k,i)*coef
              coef=coef*0.5d0*test%le
           end do
        end do
        write(10+is,*) ng, max(abs(og-test%seikai)/test%seikai,epsilon(1.d0))
     end do
  end do

  call uninit_autodiff(choose)
  call uninit_gauss_turan(maxs,maxngt,gt)
  
end program main
