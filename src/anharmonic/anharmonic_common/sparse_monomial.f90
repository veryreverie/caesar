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
    ! Constructor.
    module function new_SparseMonomial(modes) result(this) 
      type(ComplexUnivariate), intent(in) :: modes(:)
      type(SparseMonomial)                :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_SparseMonomial(this,input) 
      class(SparseMonomial), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface
    module function write_SparseMonomial(this) result(output) 
      class(SparseMonomial), intent(in) :: this
      type(String)                      :: output
    end function
  end interface
  
  interface SparseMonomial
    impure elemental module function new_SparseMonomial_String(input) &
       & result(this) 
      type(String), intent(in) :: input
      type(SparseMonomial)     :: this
    end function
  end interface
end module
