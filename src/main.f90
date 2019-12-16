program main
  use module_gauss_turan_quadrature
  implicit none

  type(GaussTuranQuadrature) :: gt(maxng)

  call init_gauss_turan(gt)

  call uninit_gauss_turan(gt)
  
end program main
