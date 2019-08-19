! ======================================================================
! A ComplexMonomial in a sparse representation to make it faster to
!    process.
! ======================================================================
module sparse_monomial_module
  use common_module
  implicit none
  
  private
  
  public :: SparseMonomial
  
  type, extends(NoDefaultConstructor) :: SparseMonomial
    type(ComplexUnivariate), allocatable :: modes(:)
  end type
  
  interface SparseMonomial
    module procedure new_SparseMonomial
  end interface
contains
function new_SparseMonomial(modes) result(this)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: modes(:)
  type(SparseMonomial)                :: this
  
  this%modes = modes
end function
end module
