! ======================================================================
! A ComplexMonomial in a sparse representation to make it faster to
!    process.
! ======================================================================
module caesar_sparse_monomial_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: SparseMonomial
  
  type, extends(Stringable) :: SparseMonomial
    type(ComplexUnivariate), allocatable :: modes(:)
  contains
    procedure, public :: read  => read_SparseMonomial
    procedure, public :: write => write_SparseMonomial
  end type
  
  interface SparseMonomial
    module procedure new_SparseMonomial
    module procedure new_SparseMonomial_String
  end interface
contains

! Constructor.
function new_SparseMonomial(modes) result(this)
  implicit none
  
  type(ComplexUnivariate), intent(in) :: modes(:)
  type(SparseMonomial)                :: this
  
  this%modes = modes
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SparseMonomial(this,input)
  implicit none
  
  class(SparseMonomial), intent(out) :: this
  type(String),          intent(in)  :: input
  
  type(ComplexUnivariate), allocatable :: modes(:)
  
  select type(this); type is(SparseMonomial)
    modes = ComplexUnivariate(tokens(input))
    this = SparseMonomial(modes)
  class default
    call err()
  end select
end subroutine

function write_SparseMonomial(this) result(output)
  implicit none
  
  class(SparseMonomial), intent(in) :: this
  type(String)                      :: output
  
  select type(this); type is(SparseMonomial)
    output = join(str(this%modes))
  class default
    call err()
  end select
end function

impure elemental function new_SparseMonomial_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(SparseMonomial)     :: this
  
  call this%read(input)
end function
end module
