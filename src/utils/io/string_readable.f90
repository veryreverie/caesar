! ======================================================================
! An abstract type, which allows extended types to be read from a String.
! ======================================================================
! Any type which extends StringReadable can be:
!    - read from String or character(*), using this=string.
!    - read as an array from an array String(:), using this(:)=string(:).
! Any type which extends StringReadable must overload %read(string).
! See example module below for how to extend this type.
module string_readable_submodule
  use io_basic_module
  use abstract_module
  
  use string_array_submodule
  implicit none
  
  private
  
  public :: StringReadable
  public :: assignment(=)
  
  type, abstract, extends(NoDefaultConstructor) :: StringReadable
  contains
    procedure(read_StringReadable), deferred :: read
  end type
  
  abstract interface
    recursive subroutine read_StringReadable(this,input)
      import String
      import StringReadable
      implicit none
      
      class(StringReadable), intent(out) :: this
      type(String),          intent(in)  :: input
    end subroutine
  end interface
  
  interface assignment(=)
    module procedure assign_StringReadable_String
    module procedure assign_StringReadable_character
  end interface
contains

! ----------------------------------------------------------------------
! Assign a StringReadable type from a String type.
! ----------------------------------------------------------------------
recursive subroutine assign_StringReadable_String(output,input)
  implicit none
  
  class(StringReadable), intent(out) :: output
  type(String),          intent(in)  :: input
  
  call output%read(input)
end subroutine

recursive subroutine assign_StringReadable_character(output,input)
  implicit none
  
  class(StringReadable), intent(out) :: output
  character(*),          intent(in)  :: input
  
  output = str(input)
end subroutine
end module

! ======================================================================
! An example module showing how to extend StringReadable.
! ======================================================================
module string_readable_example_submodule
  use io_basic_module
  
  use string_readable_submodule
  implicit none
  
  private
  
  public :: StringReadableExample
  
  type, extends(StringReadable) :: StringReadableExample
    integer :: contents
  contains
    procedure, public :: read  => read_StringReadableExample
  end type
  
  interface StringReadableExample
    module procedure new_StringReadableExample
    module procedure new_StringReadableExample_String
  end interface
contains

! Basic constructor.
function new_StringReadableExample(contents) result(this)
  implicit none
  
  integer, intent(in)         :: contents
  type(StringReadableExample) :: this
  
  this%contents = contents
end function

! The %read() routine, which is called by all StringReadable functionality.
subroutine read_StringReadableExample(this,input)
  implicit none
  
  class(StringReadableExample), intent(out) :: this
  type(String),                 intent(in)  :: input
  
  ! Select type needed to call non-polymorphic procedures,
  !    such as the StringReadableExample(input) constructor.
  select type(this); type is(StringReadableExample)
    this = StringReadableExample(int(input))
  class default
    call print_line(CODE_ERROR//': Called the StringReadableExample version &
       &of read() from a type other than StringReadableExample.')
    call err()
  end select
end subroutine

! A constructor from String, which calls %read() implicitly.
! A constructor of this format is not required,
!    but simplifies reading from file.
impure elemental function new_StringReadableExample_String(input) result(this)
  implicit none
  
  type(String), intent(in)    :: input
  type(StringReadableExample) :: this
  
  this = input
end function
end module
