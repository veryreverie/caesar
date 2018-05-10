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
module strings_readable_submodule
  use io_basic_module
  use abstract_module
  
  use string_array_submodule
  implicit none
  
  private
  
  public :: StringsReadable
  public :: assignment(=)
  
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
  
  interface assignment(=)
    module procedure assign_StringsReadable_Strings
    module procedure assign_StringsReadable_StringArray
  end interface
contains

! ----------------------------------------------------------------------
! Assign a StringsReadable type from a String type.
! ----------------------------------------------------------------------
recursive subroutine assign_StringsReadable_Strings(output,input)
  implicit none
  
  class(StringsReadable), intent(out) :: output
  type(String),           intent(in)  :: input(:)
  
  call output%read(input)
end subroutine

recursive subroutine assign_StringsReadable_StringArray(output,input)
  implicit none
  
  class(StringsReadable), intent(out) :: output
  type(StringArray),      intent(in)  :: input
  
  output = input%strings
end subroutine
end module

! ======================================================================
! An example module showing how to extend StringsReadable.
! ======================================================================
module strings_readable_example_submodule
  use io_basic_module
  
  use strings_readable_submodule
  implicit none
  
  private
  
  public :: StringsReadableExample
  
  type, extends(StringsReadable) :: StringsReadableExample
    type(String) :: line1
    type(String) :: line2
  contains
    procedure, public :: read => read_StringsReadableExample
  end type
contains

subroutine read_StringsReadableExample(this,input)
  implicit none
  
  class(StringsReadableExample), intent(out) :: this
  type(String),                  intent(in)  :: input(:)
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    read() is overloaded by any type which extends StringsReadableExample.
  select type(this); type is(StringsReadableExample)
    this%line1 = input(1)
    this%line2 = input(2)
  class default
    call print_line(CODE_ERROR//': Called the StringsReadableExample version &
       &of read() from a type other than StringsReadableExample.')
    call err()
  end select
end subroutine
end module
