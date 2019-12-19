program main
  use module_gauss_turan_quadrature
  implicit none

  type(GaussTuranQuadrature) :: gt(maxngt)

  call init_gauss_turan(gt)
  call assemble_gauss_turan
  call uninit_gauss_turan(gt)
  
end program main
