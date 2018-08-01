! ======================================================================
! An abstract type combining StringReadable and StringWriteable.
! ======================================================================
! Any type which extends Stringable has all the properties of StringReadable
!    and StringWriteable.
! Any type which extends Stringable must overload %write() and %read(string).
! See example module below for how to extend this type.
!
! N.B.: Stringable extends StringWriteable, but duplicates the functionality 
!    of StringReadable since Fortran does not support multiple inheritance.
module stringable_submodule
  use io_basic_module
  
  use string_array_submodule
  use string_readable_submodule
  use string_writeable_submodule
  implicit none
  
  private
  
  public :: Stringable
  
  type, abstract, extends(StringWriteable) :: Stringable
  contains
    procedure(read_Stringable), deferred :: read
  end type
  
  abstract interface
    recursive subroutine read_Stringable(this,input)
      import String
      import Stringable
      implicit none
      
      class(Stringable), intent(out) :: this
      type(String),      intent(in)  :: input
    end subroutine
  end interface
contains
end module

! ======================================================================
! An example module showing how to extend Stringable.
! ======================================================================
module stringable_example_submodule
  use io_basic_module
  
  use stringable_submodule
  implicit none
  
  private
  
  public :: StringableExample
  
  type, extends(Stringable) :: StringableExample
    integer :: contents
  contains
    procedure, public :: read  => read_StringableExample
    procedure, public :: write => write_StringableExample
  end type
  
  interface StringableExample
    module procedure new_StringableExample
    module procedure new_StringableExample_String
  end interface
contains

! Basic constructor.
function new_StringableExample(contents) result(this)
  implicit none
  
  integer, intent(in)     :: contents
  type(StringableExample) :: this
  
  this%contents = contents
end function

! The %read() and %write() routines,
!    as in string_readable and string_writeable.
subroutine read_StringableExample(this,input)
  implicit none
  
  class(StringableExample), intent(out) :: this
  type(String),             intent(in)  :: input
  
  select type(this); type is(StringableExample)
    this = StringableExample(int(input))
  class default
    call print_line(CODE_ERROR//': Called the StringableExample version &
       &of read() from a type other than StringableExample.')
    call err()
  end select
end subroutine

function write_StringableExample(this) result(output)
  implicit none
  
  class(StringableExample), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(StringableExample)
    output = str(this%contents)
  class default
    call print_line(CODE_ERROR//': Called the StringableExample version &
       &of write() from a type other than StringableExample.')
    call err()
  end select
end function

! The constructor from String, as in string_readable.
impure elemental function new_StringableExample_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(StringableExample)  :: this
  
  call this%read(input)
end function
end module
