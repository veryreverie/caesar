! ======================================================================
! An abstract type, which allows extended types to be read from an array
!    of strings.
! ======================================================================
! Works as StringReadable, but for types which are read from more than one
!    line.
! Any type which extends StringsReadable can be:
!    - read from String(:) or StringArray, using this=string.
!    - read as an array from StringArray(:), using this(:)=string_array(:).
! Any type which extends StringsReadable must overload %read(strings).
! See example module below or how to extend this type.
module strings_readable_module
  use io_basic_module
  use abstract_module
  
  use string_array_module
  implicit none
  
  private
  
  public :: StringsReadable
  
  type, abstract, extends(NoDefaultConstructor) :: StringsReadable
  contains
    procedure(read_StringsReadable), deferred :: read
  end type
  
  abstract interface
    recursive subroutine read_StringsReadable(this,input)
      import String
      import StringsReadable
      implicit none
      
      class(StringsReadable), intent(out) :: this
      type(String),           intent(in)  :: input(:)
    end subroutine
  end interface
contains
end module

! ======================================================================
! An example module showing how to extend StringsReadable.
! ======================================================================
module strings_readable_example_module
  use io_basic_module
  
  use string_array_module
  use strings_readable_module
  implicit none
  
  private
  
  public :: StringsReadableExample
  
  type, extends(StringsReadable) :: StringsReadableExample
    integer :: line1
    integer :: line2
  contains
    procedure, public :: read => read_StringsReadableExample
  end type
  
  interface StringsReadableExample
    module procedure new_StringsReadableExample
    module procedure new_StringsReadableExample_Strings
    module procedure new_StringsReadableExample_StringArray
  end interface
contains

! Basic constructor.
function new_StringsReadableExample(line1,line2) result(this)
  implicit none
  
  integer, intent(in)          :: line1
  integer, intent(in)          :: line2
  type(StringsReadableExample) :: this
  
  this%line1 = line1
  this%line2 = line2
end function

! The %read() subroutine.
subroutine read_StringsReadableExample(this,input)
  implicit none
  
  class(StringsReadableExample), intent(out) :: this
  type(String),                  intent(in)  :: input(:)
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    read() is overloaded by any type which extends StringsReadableExample.
  select type(this); type is(StringsReadableExample)
    this = StringsReadableExample(int(input(1)), int(input(2)))
  class default
    call print_line(CODE_ERROR//': Called the StringsReadableExample version &
       &of read() from a type other than StringsReadableExample.')
    call err()
  end select
end subroutine

! Constructor from String(:) or StringArray.
! As with string_readable, these constructor are not required but simplify
!    the read process.
! These constructors can be copied directly into user-defined classes,
!    with StringsReadableExample changed to the class name in all instances.
function new_StringsReadableExample_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)     :: input(:)
  type(StringsReadableExample) :: this
  
  call this%read(input)
end function

impure elemental function new_StringsReadableExample_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(StringsReadableExample)  :: this
  
  this = StringsReadableExample(str(input))
end function
end module
