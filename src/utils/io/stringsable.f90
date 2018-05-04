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
module stringsable_submodule
  use string_submodule
  use string_array_submodule
  use strings_readable_submodule
  use strings_writeable_submodule
  implicit none
  
  private
  
  public :: Stringsable
  public :: assignment(=)
  
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
  
  interface assignment(=)
    module procedure assign_Stringsable_Strings
    module procedure assign_Stringsable_StringArray
  end interface
contains

! ----------------------------------------------------------------------
! Assign a Stringsable type from a String type.
! ----------------------------------------------------------------------
recursive subroutine assign_Stringsable_Strings(output,input)
  implicit none
  
  class(Stringsable), intent(out) :: output
  type(String),       intent(in)  :: input(:)
  
  call output%read(input)
end subroutine

recursive subroutine assign_Stringsable_StringArray(output,input)
  implicit none
  
  class(Stringsable), intent(out) :: output
  type(StringArray),  intent(in)  :: input
  
  output = input%strings
end subroutine
end module

! ======================================================================
! An example module showing how to extend Stringsable.
! ======================================================================
module stringsable_example_submodule
  use string_submodule
  use stringsable_submodule
  use print_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: StringsableExample
  
  type, extends(Stringsable) :: StringsableExample
    type(String) :: line1
    type(String) :: line2
  contains
    procedure, public :: read => read_StringsableExample
    procedure, public :: write => write_StringsableExample
  end type
contains

subroutine read_StringsableExample(this,input)
  implicit none
  
  class(StringsableExample), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  this%line1 = input(1)
  this%line2 = input(2)
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    read() is overloaded by any type which extends StringsableExample.
  select type(this); type is(StringsableExample)
    this%line1 = input(1)
    this%line2 = input(2)
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
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    write() is overloaded by any type which extends StringsableExample.
  select type(this); type is(StringsableExample)
    output = [ this%line1, &
             & this%line2 ]
  class default
    call print_line(CODE_ERROR//': Called the StringsableExample version &
       &of write() from a type other than StringsableExample.')
    call err()
  end select
end function
end module
