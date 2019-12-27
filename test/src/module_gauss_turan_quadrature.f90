module module_gauss_turan_quadrature
  implicit none

  private
  public GaussTuranQuadrature, &
       init_gauss_turan, &
       uninit_gauss_turan 
  
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

  type(GaussTuranQuadrature),pointer :: gt(:)
  
contains
  !> allocate, associate and load gt
  subroutine init_gauss_turan(maxs,maxngt,gt_)
    integer,intent(in) :: maxs, maxngt
    type(GaussTuranQuadrature),intent(inout),target :: gt_(:)
    integer :: ig, is, iis, iin, n
    character(len=64) :: fn_gt
    
    do ig=1,maxngt
       allocate(gt_(ig)%s(0:maxs)) ! gt()%s(0)はGauss-Legendre
       do is=0,maxs
          allocate(gt_(ig)%s(is)%gz(ig)); gt_(ig)%s(is)%gz(:)=0.d0
          allocate(gt_(ig)%s(is)%we(0:2*is,ig)); gt_(ig)%s(is)%we(:,:)=0.d0
       end do
    end do
    gt=>gt_

    write(fn_gt, '(a,i3.3,a,i3.3,a)') "gt_s",maxs,"_n",maxngt,".conf"    
    write(*,*) fn_gt
    open(1,file=fn_gt)
    do iis=0,maxs
       read(1,*) is
       if(is.ne.iis) stop "aho"
       do iin=1,maxngt
          read(1,*) n
          if(iin.ne.n) stop "aho"
          do ig=1,n
             read(1,*) gt(n)%s(is)%gz(ig), gt(n)%s(is)%we(0:2*is,ig)
          end do
       end do
    end do
    close(1)

  end subroutine init_gauss_turan

  !> deallocate and nullify gt  
  subroutine uninit_gauss_turan(maxs,maxngt,gt_)
    integer,intent(in) :: maxs, maxngt
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
  
end module module_gauss_turan_quadrature
