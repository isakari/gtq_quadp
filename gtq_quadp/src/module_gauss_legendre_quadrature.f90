module module_gauss_legendre_quadrature
  use module_qlapack
  implicit none

  private
  public assemble_gauss_legendre

contains
  !> ガウス積分の分点と重みを計算する
  !! \param ng 分点の数
  !! \param gz 分点の座標
  !! \param we 重み
  subroutine assemble_gauss_legendre(n,gz,we)
    integer,intent(in) :: n
    real(16),intent(out) :: gz(n), we(n)

    ! 対称三重対角行列の対角成分
    real(16) :: td(n)

    ! 対称三重対角行列の非対角成分
    real(16) :: tl(n-1)

    ! 固有ベクトル
    complex(16) :: q(n,n)
    
    integer :: i, info
    real(16) :: work(2*n-2)
    
    td(:)=0.q0
    do i=1,n-1
       tl(i)=dble(i)/sqrt((2.q0*i-1.0q0)*(2.q0*i+1.q0))
    end do

    q(:,:)=cmplx(0.q0,0.q0,kind(1.q0))
    do i=1,n
       q(i,i)=cmplx(1.q0,0.q0,kind(1.q0))
    end do

    call qzsteqr("V",n,td,tl,q,n,work,info)

    do i=1,n
       gz(i)=td(i)
       we(i)=2.q0*(real(q(1,i))**2+aimag(q(1,i))**2)
    end do
    
  end subroutine assemble_gauss_legendre
  
end module module_gauss_legendre_quadrature
  
