! ======================================================================
! A vector along a single complex mode.
! ======================================================================
module complex_single_mode_vector_submodule
  use utils_module
  implicit none
  
  private
  
  public :: ComplexSingleModeVector
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
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
  
  interface operator(*)
    module procedure multiply_real_ComplexSingleModeVector
    module procedure multiply_ComplexSingleModeVector_real
    module procedure multiply_complex_ComplexSingleModeVector
    module procedure multiply_ComplexSingleModeVector_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexSingleModeVector_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexSingleModeVector_ComplexSingleModeVector
  end interface
  
  interface operator(-)
    module procedure negative_ComplexSingleModeVector
    module procedure subtract_ComplexSingleModeVector_ComplexSingleModeVector
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

! Arithmetic.
impure elemental function multiply_real_ComplexSingleModeVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),                      intent(in) :: this
  type(ComplexSingleModeVector), intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector( id        = that%id, &
                                  & magnitude = this*that%magnitude)
end function

impure elemental function multiply_ComplexSingleModeVector_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  real(dp),                      intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector( id        = this%id, &
                                  & magnitude = this%magnitude*that)
end function

impure elemental function multiply_complex_ComplexSingleModeVector(this,that) &
   & result(output)
  implicit none
  
  complex(dp),                   intent(in) :: this
  type(ComplexSingleModeVector), intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector( id        = that%id, &
                                  & magnitude = this*that%magnitude)
end function

impure elemental function multiply_ComplexSingleModeVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector( id        = this%id, &
                                  & magnitude = this%magnitude*that)
end function

impure elemental function divide_ComplexSingleModeVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector( id        = this%id, &
                                  & magnitude = this%magnitude/that)
end function

impure elemental function add_ComplexSingleModeVector_ComplexSingleModeVector(&
   & this,that) result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  type(ComplexSingleModeVector), intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleModeVector( id        = this%id, &
                                  & magnitude = this%magnitude+that%magnitude)
end function

impure elemental function negative_ComplexSingleModeVector(this) result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  type(ComplexSingleModeVector)             :: output
  
  output = ComplexSingleModeVector(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function                                                &
   & subtract_ComplexSingleModeVector_ComplexSingleModeVector(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleModeVector), intent(in) :: this
  type(ComplexSingleModeVector), intent(in) :: that
  type(ComplexSingleModeVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleModeVector( id        = this%id, &
                                  & magnitude = this%magnitude-that%magnitude)
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
