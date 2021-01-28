! ======================================================================
! An abstract type combining StringsReadable and StringsWriteable.
! ======================================================================
! Any type which extends Stringsable has all the properties of StringsReadable
!    and StringsWriteable.
! Any type which extends Stringsable must overload %write() and %read(strings).
! See example module below or how to extend this type.
!
! N.B.: Stringsable extends StringsWriteable, but duplicates the functionality 
!    of StringsReadable since Fortran does not support multiple inheritance.
module caesar_stringsable_module
  use caesar_io_basic_module
  
  use caesar_string_array_module
  use caesar_strings_readable_module
  use caesar_strings_writeable_module
  implicit none
  
  private
  
  public :: Stringsable
  
  type, abstract, extends(StringsWriteable) :: Stringsable
  contains
    procedure(read_Stringsable), deferred :: read
  end type
  
  abstract interface
    recursive subroutine read_Stringsable(this,input)
      import String
      import Stringsable
      implicit none
      
      class(Stringsable), intent(out) :: this
      type(String),       intent(in)  :: input(:)
    end subroutine
  end interface
contains
end module

! ======================================================================
! An example module showing how to extend Stringsable.
! ======================================================================
module caesar_stringsable_example_module
  use caesar_io_basic_module
  
  use caesar_string_array_module
  use caesar_stringsable_module
  implicit none
  
  private
  
  public :: StringsableExample
  
  type, extends(Stringsable) :: StringsableExample
    integer :: line1
    integer :: line2
  contains
    procedure, public :: read => read_StringsableExample
    procedure, public :: write => write_StringsableExample
  end type
  
  interface StringsableExample
    module procedure new_StringsableExample
    module procedure new_StringsableExample_Strings
    module procedure new_StringsableExample_StringArray
  end interface
contains

! Basic constructor.
function new_StringsableExample(line1,line2) result(this)
  implicit none
  
  integer, intent(in)      :: line1
  integer, intent(in)      :: line2
  type(StringsableExample) :: this
  
  this%line1 = line1
  this%line2 = line2
end function

! The %read() and %write() routines,
!    as in strings_readable and strings_writeable.
subroutine read_StringsableExample(this,input)
  implicit none
  
  class(StringsableExample), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  select type(this); type is(StringsableExample)
    this = StringsableExample(int(input(1)), int(input(2)))
  class default
    call print_line(CODE_ERROR//': Called the StringsableExample version &
       &of read() from a type other than StringsableExample.')
    call err()
  end select
end subroutine

function write_StringsableExample(this) result(output)
  implicit none
  
  class(StringsableExample), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(StringsableExample)
    output = [ str(this%line1), &
             & str(this%line2) ]
  class default
    call print_line(CODE_ERROR//': Called the StringsableExample version &
       &of write() from a type other than StringsableExample.')
    call err()
  end select
end function

! Constructor from String(:) or StringArray.
! As with string_readable, these constructor are not required but simplify
!    the read process.
! These constructors can be copied directly into user-defined classes,
!    with StringsableExample changed to the class name in all instances.
function new_StringsableExample_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StringsableExample) :: this
  
  call this%read(input)
end function

impure elemental function new_StringsableExample_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StringsableExample)      :: this
  
  this = StringsableExample(str(input))
end function
end module
