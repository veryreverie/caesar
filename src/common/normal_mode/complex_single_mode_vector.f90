! ======================================================================
! A vector along a single complex mode.
! ======================================================================
module complex_single_mode_vector_submodule
  use utils_module
  implicit none
  
  private
  
  public :: ComplexSingleModeVector
  
  type, extends(Stringable) :: ComplexSingleModeVector
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the vector along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleModeVector
    procedure, public :: write => write_ComplexSingleModeVector
  end type
  
  interface ComplexSingleModeVector
    module procedure new_ComplexSingleModeVector
    module procedure new_ComplexSingleModeVector_String
  end interface
contains

! Constructor.
function new_ComplexSingleModeVector(id,magnitude) result(this)
  implicit none
  
  integer,     intent(in)       :: id
  complex(dp), intent(in)       :: magnitude
  type(ComplexSingleModeVector) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

! I/O.
subroutine read_ComplexSingleModeVector(this,input)
  implicit none
  
  class(ComplexSingleModeVector), intent(out) :: this
  type(String),                   intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleModeVector)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse complex single mode &
         &vector from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1+1.2i then split_string = ["u3","=","2.1+1.2i"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    magnitude = cmplx(split_string(3))
    
    this = ComplexSingleModeVector(id,magnitude)
  end select
end subroutine

function write_ComplexSingleModeVector(this) result(output)
  implicit none
  
  class(ComplexSingleModeVector), intent(in) :: this
  type(String)                               :: output
  
  select type(this); type is(ComplexSingleModeVector)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_ComplexSingleModeVector_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)      :: input
  type(ComplexSingleModeVector) :: this
  
  this = input
end function
end module
