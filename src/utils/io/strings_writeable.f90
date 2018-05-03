! ======================================================================
! An abstract type, which allows extended types to be written to an array
!    of strings.
! ======================================================================
! Works as StringWriteable, but for types which require more than one line
!    to print.
! Any type which extends StringsWriteable can be:
!    - converted to String(:) using str(this).
!    - printed to stdout, using print_lines(this).
!    - printed to file, using file%print_lines(this).
! An array of StringWriteables can be printed across multiple lines using
!    print_lines(this(:)) or file%print_lines(this(:)),
!    with the optional keyword separating_line to separate elements.
! Any type which extends StringsWriteable must overload %write().
! See example module below for how to extend this type.
module strings_writeable_submodule
  use string_submodule
  use string_array_submodule
  use error_submodule
  use print_submodule
  implicit none
  
  private
  
  public :: StringsWriteable
  public :: str
  public :: print_lines
  
  type, abstract :: StringsWriteable
  contains
    procedure(write_StringsWriteable), deferred :: write
  end type
  
  abstract interface
    recursive function write_StringsWriteable(this) result(output)
      import String
      import StringsWriteable
      implicit none
      
      class(StringsWriteable), intent(in) :: this
      type(String), allocatable           :: output(:)
    end function
  end interface
  
  interface str
    module procedure str_StringsWriteable
  end interface
  
  interface print_lines
    module procedure print_lines_StringsWriteable
    module procedure print_lines_StringsWriteables_character
    module procedure print_lines_StringsWriteables_String
  end interface
contains

! ----------------------------------------------------------------------
! The str() function, which converts to string types.
! ----------------------------------------------------------------------
recursive function str_StringsWriteable(this) result(output)
  implicit none
  
  class(StringsWriteable), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  output = this%write()
end function

! ----------------------------------------------------------------------
! Provides print_lines() for types which extend StringsWriteable.
! ----------------------------------------------------------------------
subroutine print_lines_StringsWriteable(this,indent)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this
  integer,                 intent(in), optional :: indent
  
  call print_lines(str(this), indent)
end subroutine

subroutine print_lines_StringsWriteables_character(this,indent,separating_line)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  integer,                 intent(in), optional :: indent
  character(*),            intent(in), optional :: separating_line
  
  integer :: i
  
  do i=1,size(this)
    call print_lines(this(i), indent)
    if (present(separating_line)) then
      call print_line(separating_line, indent)
    endif
  enddo
end subroutine

subroutine print_lines_StringsWriteables_String(this,indent,separating_line)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  integer,                 intent(in), optional :: indent
  type(String),            intent(in)           :: separating_line
  
  call print_lines(this, indent, char(separating_line))
end subroutine
end module

! ======================================================================
! An example module showing how to extend the StringsWriteable type.
! ======================================================================
module StringsWriteable_example_submodule
  use string_submodule
  use strings_writeable_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: StringsWriteableExample
  
  type, extends(StringsWriteable) :: StringsWriteableExample
    type(String) :: line1
    type(String) :: line2
  contains
    procedure :: write => write_StringsWriteableExample
  end type
contains

function write_StringsWriteableExample(this) result(output)
  implicit none
  
  class(StringsWriteableExample), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    write() is overloaded by any type which extends StringsWriteableExample.
  select type(this); type is(StringsWriteableExample)
    output = [ this%line1, &
             & this%line2 ]
  class default
    call print_line(CODE_ERROR//': Called the StringsWriteableExample version &
       &of write() from a type other than StringsWriteableExample.')
    call err()
  end select
end function
end module
