module module_gauss_turan_quadrature
  implicit none

  private
  public GaussTuranQuadrature, &
       maxs, maxng, &
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
  
  integer,parameter :: maxs=2 !< maximum s
  integer,parameter :: maxng=3 !< maximum number of nodes

  type(GaussTuranQuadrature),pointer :: gt(:)

contains
  !> allocate and associate gt
  subroutine init_gauss_turan(gt_)
    type(GaussTuranQuadrature),intent(inout),target :: gt_(:)
    integer :: ig, is

    do ig=1,maxng
       allocate(gt_(ig)%s(0:maxs))
       do is=0,maxs
          allocate(gt_(ig)%s(is)%gz(ig))
          allocate(gt_(ig)%s(is)%we(0:2*is,ig))
       end do
    end do

    gt=>gt_
    
  end subroutine init_gauss_turan

  !> deallocate and dissociate gt  
  subroutine uninit_gauss_turan(gt_)
    type(GaussTuranQuadrature),intent(inout),target :: gt_(:)
    integer :: ig, is

    do ig=1,maxng
       do is=0,maxs
          if(allocated(gt_(ig)%s(is)%gz)) deallocate(gt_(ig)%s(is)%gz)
          if(allocated(gt_(ig)%s(is)%we)) deallocate(gt_(ig)%s(is)%we)
       end do
       if(allocated(gt_(ig)%s)) deallocate(gt_(ig)%s)
    end do

    if(associated(gt)) nullify(gt)
    
  end subroutine uninit_gauss_turan
  
end module module_gauss_turan_quadrature
