module module_gauss_turan_quadrature
  use module_math
  use module_gauss_legendre_quadrature
  implicit none

  private
  public GaussTuranQuadrature, &
       maxs, maxngt, &
       init_gauss_turan, &
       uninit_gauss_turan, &
       assemble_gauss_turan

  type GaussTuranQuadrature2
     integer :: ng !< the num of nodes
     real(8),allocatable :: gz(:) !< nodes
     real(8),allocatable :: we(:,:) !< weights
  end type GaussTuranQuadrature2

  !> type(GaussTuranQuadrature) gt(n)%s(s) contains nodes and weights for Gauss-Turan-type quadrature (GTQ)
  !> which uses function value and derivatives (up to 2s-th order) at n nodes.
  type GaussTuranQuadrature
     type(GaussTuranQuadrature2),allocatable :: s(:)
  end type GaussTuranQuadrature

  integer,parameter :: maxs=2 !< maximum s
  integer,parameter :: maxngt=3 !< maximum number of nodes
  
  !> type(GaussLegendreQuadrature) gl(n) contains nodes and weights for Gauss-Legendre quadrature (GLQ)
  !> which uses function value at n nodes.
  type GaussLegendreQuadrature
     integer :: ng !< the num of nodes
     real(8),allocatable :: gz(:) !< nodes
     real(8),allocatable :: we(:) !< weights
  end type GaussLegendreQuadrature

  integer,parameter :: maxngl=(maxs+1)*maxngt
  type(GaussLegendreQuadrature) :: gl(maxngl)

  ! pointer
  type(GaussTuranQuadrature),pointer :: gt(:)

contains
  !> allocate and associate gt
  subroutine init_gauss_turan(gt_)
    type(GaussTuranQuadrature),intent(inout),target :: gt_(:)
    integer :: ig, is

    do ig=1,maxngt
       allocate(gt_(ig)%s(0:maxs)) ! gt()%s(0)はGauss-Legendre
       do is=0,maxs
          allocate(gt_(ig)%s(is)%gz(ig)); gt_(ig)%s(is)%gz(:)=0.d0
          allocate(gt_(ig)%s(is)%we(0:2*is,ig)); gt_(ig)%s(is)%we(:,:)=0.d0
       end do
    end do

    gt=>gt_
    
  end subroutine init_gauss_turan

  !> deallocate and dissociate gt  
  subroutine uninit_gauss_turan(gt_)
    type(GaussTuranQuadrature),intent(inout),target :: gt_(:)

    integer :: ig, is

    do ig=1,maxngt
       do is=0,maxs
          if(allocated(gt_(ig)%s(is)%gz)) deallocate(gt_(ig)%s(is)%gz)
          if(allocated(gt_(ig)%s(is)%we)) deallocate(gt_(ig)%s(is)%we)
       end do
       if(allocated(gt_(ig)%s)) deallocate(gt_(ig)%s)
    end do

    if(associated(gt)) nullify(gt)
    
  end subroutine uninit_gauss_turan

  !> Kroneckerのデルタ
  real(8) function delta(i,j) result(out)
    integer,intent(in) :: i, j
    if(i.eq.j) then
       out=1.d0
    else
       out=0.d0
    end if
  end function delta

  !> Turan多項式
  recursive real(8) function poly(n,i,bet,x) result(out)
    integer,intent(in) :: n
    integer,intent(in) :: i
    real(8),intent(in) :: bet(n-1)
    real(8),intent(in) :: x
    if(i.eq.-1) then
       out=0.d0
    elseif(i.eq.0) then
       out=1.d0
    elseif(i.eq.1)then
       out=x
    else
       out=x*poly(n,i-1,bet,x)-bet(i-1)*poly(n,i-2,bet,x)
    end if
  end function poly

  !> Turan多項式のbeta微分
  recursive real(8) function dpoly(n,i,j,bet,x) result(out)
    integer,intent(in) :: n
    integer,intent(in) :: i, j
    real(8),intent(in) :: bet(n-1)
    real(8),intent(in) :: x
    if(i.le.j) then
       out=0.d0
    elseif(i.eq.j+1)then
       out=-poly(n,j-1,bet,x)
    else
       out=x*dpoly(n,i-1,j,bet,x)-bet(i-1)*dpoly(n,i-2,j,bet,x)
    end if
  end function dpoly
  
  !> compute nodes and weights for GTQ
  subroutine assemble_gauss_turan
    integer :: i, j, k, l, is, n, ngl, ig, ipiv(maxngt-1), info, itmp
    real(8) :: bet(maxngt-1), w(maxngt-1), jac(maxngt-1,maxngt-1), f(maxngt-1)
    real(8) :: p(0:maxngt), q(0:maxngt,0:maxngt), b(0:maxngt,0:maxngt)
    real(8) :: td(maxngt), tl(maxngt-1), work(2*maxngt-2)
    complex(8) :: ev(maxngt,maxngt)

    real(8) :: ahat(2*maxs+1,2*maxs+1), u(2*maxs), muhat(0:2*maxs), tmp, bb(2*maxs+1)

!!$    real(8) :: x
    
    ! GLを準備
    do i=1,maxngl
       gl(i)%ng=i
       allocate(gl(i)%gz(i),gl(i)%we(i))
       call assemble_gauss_legendre(i,gl(i)%gz,gl(i)%we)
    end do

    ! 分点数をset
    do is=0,maxs
       do i=1,maxngt
          gt(i)%s(is)%ng=i
       end do
    end do
    
    ! s=0解微分を使うGauss-TuranはGauss-Legendreに他ならない
    do i=1,maxngt
       gt(i)%s(0)%gz(:)=gl(i)%gz(:)
       gt(i)%s(0)%we(0,:)=gl(i)%we(:)
    end do
    
    do is=1,maxs !is=(1,..,maxs)階微分を使うGauss-Turan多項式の零点を探す
       
       ! Turan多項式の満たす漸化式の係数(beta)をnewton法で求める
       ! P(n;-1)(x)=0
       ! P(n;0)(x)=1
       ! P(n;1)(x)=x
       ! P(n;k+1)(x)=xP(n;k)(x)-beta(k)P(n;k-1)(x)

       ! n=2のときのNewton法の初期値
       bet(1)=(2.d0*is+1.d0)/(2.d0*is+3.d0)
       bet(3:)=0.d0

       do n=2,maxngt ! n=2次公式から順番に
          ngl=(is+1)*n ! ngl点のGL公式でfの値を計算し
          if(n.gt.2) bet(n-1)=bet(n-2)
          do itmp=1,100! Newton法で
             f(:)=0.d0 ! f(beta)=0となるbetaを探す
             jac(:,:)=0.d0 ! <=Newton法で使うJacobi行列
             ! fとそのJacobianを計算
             do ig=1,ngl
                do i=0,n
                   p(i)=poly(n,i,bet,gl(ngl)%gz(ig))
                   do j=0,n
                      b(i,j)=dpoly(n,i,j,bet,gl(ngl)%gz(ig))
                   end do
                end do
                do i=0,n
                   do j=0,n
                      q(i,j)=p(i)*(b(i,j)*p(n)+is*b(n,j)*p(i))
                   end do
                end do
                do i=1,n-1
                   f(i)=f(i)+(bet(i)*p(i-1)**2-p(i)**2)*p(n)**(2*is)*gl(ngl)%we(ig)
                   do j=1,n-1
                      jac(i,j)=jac(i,j)+2.d0*(bet(i)*q(i-1,j)-q(i,j)&
                           +0.5d0*delta(i,j)*p(i-1)**2*p(n))*p(n)**(2*is-1)*gl(ngl)%we(ig)
                   end do
                end do
             end do
             ! Newton法の更新量を求める
             w=f
             call dgesv(n-1,1,jac(1:n-1,1:n-1),n-1,ipiv(1:n-1),w(1:n-1),n-1,info)
             ! write(100+n,*) bet(1:n-1)
             if(sqrt(dot_product(w(1:n-1),w(1:n-1))).le.epsilon(1.d0)) exit
             bet(:)=bet(:)-w(:) ! 更新
          end do
          ! write(*,'(2i3,100f24.18)') is,n,bet(1:n-1)

          td(:)=0.d0
          do i=1,n-1
             tl(i)=sqrt(bet(i))
          end do

          ev(:,:)=cmplx(0.d0,0.d0,kind(1.d0))
          do i=1,n
             ev(i,i)=cmplx(1.d0,0.d0,kind(1.d0))
          end do

          call zsteqr("V",n,td(1:n),tl(1:n-1),ev(1:n,1:n),n,work(1:n),info)
          ! write(*,*) "info", info
          
          do i=1,n
             gt(n)%s(is)%gz(i)=td(i)
          end do
          
       end do

       ! n=1の時は分点は零
       gt(1)%s(is)%gz(1)=0.d0

    end do

    do is=1,maxs ! 2*is階微分を使うGT公式の
       do n=1,maxngt ! n次の公式の
          do j=1,n ! j番分点の重みを計算する

             ! uを計算
             u(:)=0.d0
             do i=1,2*is
                do l=1,n
                   if(l.ne.j) then
                      u(i)=u(i)+(gt(n)%s(is)%gz(l)-gt(n)%s(is)%gz(j))**(-i)
                   end if
                end do
             end do

             ! ahatを計算
             ahat(:,:)=0.d0
             do k=1,2*is+1
                ahat(k,k)=1.d0
             end do
             do l=1,2*is
                do k=1,2*is+1-l
                   do i=1,l
                      ahat(k,k+l)=ahat(k,k+l)-(2.d0*is+1.d0)/dble(l)*u(i)*ahat(i,l)
                   end do
                end do
             end do

             ! muhatを計算
             ngl=(is+1)*n
             muhat(:)=0.d0
             do k=0,2*is
                do ig=1,ngl
                   tmp=1.d0
                   do l=1,n
                      if(l.ne.j) then
                         tmp=tmp*((gl(ngl)%gz(ig)-gt(n)%s(is)%gz(l))&
                              /(gt(n)%s(is)%gz(j)-gt(n)%s(is)%gz(l)))**(2*is+1)
                      end if
                   end do
                   muhat(k)=muhat(k)+(gl(ngl)%gz(ig)-gt(n)%s(is)%gz(j))**k*tmp*gl(ngl)%we(ig)
                end do
             end do

             ! bを計算
             bb(:)=0.d0
             bb(2*is+1)=muhat(2*is)
             do k=2*is,1,-1
                bb(k)=muhat(k-1)
                do l=k+1,2*is+1
                   bb(k)=bb(k)-ahat(k,l)*bb(l)
                end do
             end do

             ! bb-->weight
             do k=1,2*is+1
                gt(n)%s(is)%we(k-1,j)=bb(k)/factorial(k-1)
             end do
             
          end do
       end do
    end do

    ! 分点・重みをcheck
    do is=0,maxs
       do n=1,maxngt
          write(*,*) "#", n,"-points GTQ with up to",2*is,"-th derivatives"
          write(*,*) "# i, i-th node (x_i), weight for f^(s)(x_i), s=0,...",2*is
          do i=1,n
             write(*,*) i,gt(n)%s(is)%gz(i), gt(n)%s(is)%we(0:2*is,i)
          end do
       end do
    end do

!!$    write(*,*) "s",sin(1.d0)-sin(-1.d0)
!!$
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(0)%ng
!!$       x=gt(10)%s(1)%gz(ig)
!!$       tmp=tmp+cos(x)*gt(10)%s(0)%we(0,ig)
!!$    end do
!!$    write(*,*) "0",tmp
!!$    
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(1)%ng
!!$       x=gt(10)%s(1)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(1)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(1)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(1)%we(2,ig))
!!$    end do
!!$    write(*,*) "1",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(2)%ng
!!$       x=gt(10)%s(2)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(2)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(2)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(2)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(2)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(2)%we(4,ig))
!!$    end do
!!$    write(*,*) "2",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(3)%ng
!!$       x=gt(10)%s(3)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(3)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(3)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(3)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(3)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(3)%we(4,ig)&
!!$               -sin(x)*gt(10)%s(3)%we(5,ig)&
!!$               -cos(x)*gt(10)%s(3)%we(6,ig))
!!$    end do
!!$    write(*,*) "3",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(4)%ng
!!$       x=gt(10)%s(4)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(4)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(4)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(4)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(4)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(4)%we(4,ig)&
!!$               -sin(x)*gt(10)%s(4)%we(5,ig)&
!!$               -cos(x)*gt(10)%s(4)%we(6,ig)&
!!$               +sin(x)*gt(10)%s(4)%we(7,ig)&
!!$               +cos(x)*gt(10)%s(4)%we(8,ig))
!!$    end do
!!$    write(*,*) "4",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(5)%ng
!!$       x=gt(10)%s(5)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(5)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(5)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(5)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(5)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(5)%we(4,ig)&
!!$               -sin(x)*gt(10)%s(5)%we(5,ig)&
!!$               -cos(x)*gt(10)%s(5)%we(6,ig)&
!!$               +sin(x)*gt(10)%s(5)%we(7,ig)&
!!$               +cos(x)*gt(10)%s(5)%we(8,ig)&
!!$               -sin(x)*gt(10)%s(5)%we(9,ig)&
!!$               -cos(x)*gt(10)%s(5)%we(10,ig))
!!$    end do
!!$    write(*,*) "5",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(6)%ng
!!$       x=gt(10)%s(6)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(6)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(6)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(6)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(6)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(6)%we(4,ig)&
!!$               -sin(x)*gt(10)%s(6)%we(5,ig)&
!!$               -cos(x)*gt(10)%s(6)%we(6,ig)&
!!$               +sin(x)*gt(10)%s(6)%we(7,ig)&
!!$               +cos(x)*gt(10)%s(6)%we(8,ig)&
!!$               -sin(x)*gt(10)%s(6)%we(9,ig)&
!!$               -cos(x)*gt(10)%s(6)%we(10,ig)&
!!$               +sin(x)*gt(10)%s(6)%we(11,ig)&
!!$               +cos(x)*gt(10)%s(6)%we(12,ig))
!!$    end do
!!$    write(*,*) "6",tmp
!!$    tmp=0.d0
!!$    do ig=1,gt(10)%s(7)%ng
!!$       x=gt(10)%s(7)%gz(ig)
!!$       tmp=tmp+(cos(x)*gt(10)%s(7)%we(0,ig)&
!!$               -sin(x)*gt(10)%s(7)%we(1,ig)&
!!$               -cos(x)*gt(10)%s(7)%we(2,ig)&
!!$               +sin(x)*gt(10)%s(7)%we(3,ig)&
!!$               +cos(x)*gt(10)%s(7)%we(4,ig)&
!!$               -sin(x)*gt(10)%s(7)%we(5,ig)&
!!$               -cos(x)*gt(10)%s(7)%we(6,ig)&
!!$               +sin(x)*gt(10)%s(7)%we(7,ig)&
!!$               +cos(x)*gt(10)%s(7)%we(8,ig)&
!!$               -sin(x)*gt(10)%s(7)%we(9,ig)&
!!$               -cos(x)*gt(10)%s(7)%we(10,ig)&
!!$               +sin(x)*gt(10)%s(7)%we(11,ig)&
!!$               +cos(x)*gt(10)%s(7)%we(12,ig))
!!$    end do
!!$    write(*,*) "7",tmp
    
    ! GLを捨てる
    do i=1,maxngl
       if(allocated(gl(i)%gz)) deallocate(gl(i)%gz)
       if(allocated(gl(i)%we)) deallocate(gl(i)%we)
    end do
    
  end subroutine assemble_gauss_turan
  
end module module_gauss_turan_quadrature
