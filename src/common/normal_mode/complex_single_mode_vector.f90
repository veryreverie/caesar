! ======================================================================
! A vector along a single complex mode.
! ======================================================================
module complex_single_mode_vector_submodule
  use utils_module
  
  use complex_mode_submodule
  implicit none
  
  private
  
  public :: ComplexSingleVector
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  
  type, extends(Stringable) :: ComplexSingleVector
    ! The id of the mode.
    integer :: id
    
    ! The magnitude of the vector along the mode.
    complex(dp) :: magnitude
  contains
    procedure, public :: read  => read_ComplexSingleVector
    procedure, public :: write => write_ComplexSingleVector
  end type
  
  interface ComplexSingleVector
    module procedure new_ComplexSingleVector
    module procedure new_ComplexSingleVector_ComplexMode
    module procedure new_ComplexSingleVector_String
  end interface
  
  interface operator(*)
    module procedure multiply_real_ComplexSingleVector
    module procedure multiply_ComplexSingleVector_real
    module procedure multiply_complex_ComplexSingleVector
    module procedure multiply_ComplexSingleVector_complex
  end interface
  
  interface operator(/)
    module procedure divide_ComplexSingleVector_complex
  end interface
  
  interface operator(+)
    module procedure add_ComplexSingleVector_ComplexSingleVector
  end interface
  
  interface operator(-)
    module procedure negative_ComplexSingleVector
    module procedure subtract_ComplexSingleVector_ComplexSingleVector
  end interface
contains

! Constructors.
impure elemental function new_ComplexSingleVector(id,magnitude) &
   & result(this)
  implicit none
  
  integer,     intent(in)       :: id
  complex(dp), intent(in)       :: magnitude
  type(ComplexSingleVector) :: this
  
  this%id        = id
  this%magnitude = magnitude
end function

impure elemental function new_ComplexSingleVector_ComplexMode(mode, &
   & magnitude) result(this)
  implicit none
  
  type(ComplexMode), intent(in) :: mode
  complex(dp),       intent(in) :: magnitude
  type(ComplexSingleVector) :: this
  
  this = ComplexSingleVector(id=mode%id, magnitude=magnitude)
end function

! Arithmetic.
impure elemental function multiply_real_ComplexSingleVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),                      intent(in) :: this
  type(ComplexSingleVector), intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector( id        = that%id, &
                                  & magnitude = this*that%magnitude)
end function

impure elemental function multiply_ComplexSingleVector_real(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  real(dp),                      intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector( id        = this%id, &
                                  & magnitude = this%magnitude*that)
end function

impure elemental function multiply_complex_ComplexSingleVector(this,that) &
   & result(output)
  implicit none
  
  complex(dp),                   intent(in) :: this
  type(ComplexSingleVector), intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector( id        = that%id, &
                                  & magnitude = this*that%magnitude)
end function

impure elemental function multiply_ComplexSingleVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector( id        = this%id, &
                                  & magnitude = this%magnitude*that)
end function

impure elemental function divide_ComplexSingleVector_complex(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  complex(dp),                   intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector( id        = this%id, &
                                  & magnitude = this%magnitude/that)
end function

impure elemental function add_ComplexSingleVector_ComplexSingleVector(&
   & this,that) result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  type(ComplexSingleVector), intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to add vectors along different modes.')
    call err()
  endif
  
  output = ComplexSingleVector( id        = this%id, &
                                  & magnitude = this%magnitude+that%magnitude)
end function

impure elemental function negative_ComplexSingleVector(this) result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  type(ComplexSingleVector)             :: output
  
  output = ComplexSingleVector(id=this%id, magnitude=-this%magnitude)
end function

impure elemental function                                                &
   & subtract_ComplexSingleVector_ComplexSingleVector(this,that) &
   & result(output)
  implicit none
  
  type(ComplexSingleVector), intent(in) :: this
  type(ComplexSingleVector), intent(in) :: that
  type(ComplexSingleVector)             :: output
  
  if (this%id/=that%id) then
    call print_line(ERROR//': Trying to subtract vectors along different &
       &modes.')
    call err()
  endif
  
  output = ComplexSingleVector( id        = this%id, &
                                  & magnitude = this%magnitude-that%magnitude)
end function

! I/O.
subroutine read_ComplexSingleVector(this,input)
  implicit none
  
  class(ComplexSingleVector), intent(out) :: this
  type(String),                   intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: magnitude
  
  select type(this); type is(ComplexSingleVector)
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
    
    this = ComplexSingleVector(id,magnitude)
  end select
end subroutine

function write_ComplexSingleVector(this) result(output)
  implicit none
  
  class(ComplexSingleVector), intent(in) :: this
  type(String)                               :: output
  
  select type(this); type is(ComplexSingleVector)
    output = 'u'//this%id//' = '//this%magnitude
  end select
end function

impure elemental function new_ComplexSingleVector_String(input) &
   & result(this)
  implicit none
  
  type(String), intent(in)      :: input
  type(ComplexSingleVector) :: this
  
  this = input
end function
end module
