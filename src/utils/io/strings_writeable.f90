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
  use io_basic_module
  use abstract_module
  
  use string_array_submodule
  implicit none
  
  private
  
  public :: StringsWriteable
  public :: str
  public :: print_lines
  
  type, abstract, extends(NoDefaultConstructor) :: StringsWriteable
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
    module procedure str_StringsWriteables_String
    module procedure str_StringsWriteables_character
  end interface
  
  interface print_lines
    module procedure print_lines_StringsWriteable
    module procedure print_lines_StringsWriteables_String
    module procedure print_lines_StringsWriteables_character
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

recursive function str_StringsWriteables_String(this,separating_line) &
   & result(output)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  type(String),            intent(in), optional :: separating_line
  type(String), allocatable                     :: output(:)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,ialloc
  
  allocate(sections(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    sections(i) = StringArray(str(this(i)))
  enddo
  
  output = str(join(sections,separating_line))
end function

recursive function str_StringsWriteables_character(this,separating_line) &
   & result(output)
  implicit none
  
  class(StringsWriteable), intent(in) :: this(:)
  character(*),            intent(in) :: separating_line
  type(String), allocatable           :: output(:)
  
  output = str(this,str(separating_line))
end function

! ----------------------------------------------------------------------
! Provides print_lines() for types which extend StringsWriteable.
! ----------------------------------------------------------------------
subroutine print_lines_StringsWriteable(this)
  implicit none
  
  class(StringsWriteable), intent(in) :: this
  
  call print_lines(str(this))
end subroutine

subroutine print_lines_StringsWriteables_String(this,separating_line)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  type(String),            intent(in), optional :: separating_line
  
  call print_lines(str(this,separating_line))
end subroutine

subroutine print_lines_StringsWriteables_character(this,separating_line)
  implicit none
  
  class(StringsWriteable), intent(in) :: this(:)
  character(*),            intent(in) :: separating_line
  
  call print_lines(str(this,separating_line))
end subroutine
end module

! ======================================================================
! An example module showing how to extend the StringsWriteable type.
! ======================================================================
module StringsWriteable_example_submodule
  use io_basic_module
  
  use strings_writeable_submodule
  implicit none
  
  private
  
  public :: StringsWriteableExample
  
  type, extends(StringsWriteable) :: StringsWriteableExample
    integer :: line1
    integer :: line2
  contains
    procedure :: write => write_StringsWriteableExample
  end type
contains

! The %write() function.
function write_StringsWriteableExample(this) result(output)
  implicit none
  
  class(StringsWriteableExample), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    write() is overloaded by any type which extends StringsWriteableExample.
  select type(this); type is(StringsWriteableExample)
    output = [ str(this%line1), &
             & str(this%line2) ]
  class default
    call print_line(CODE_ERROR//': Called the StringsWriteableExample version &
       &of write() from a type other than StringsWriteableExample.')
    call err()
  end select
end function
end module
