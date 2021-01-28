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
module caesar_strings_writeable_module
  use caesar_io_basic_module
  use caesar_abstract_module
  
  use caesar_string_array_module
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
recursive function str_StringsWriteable(this,settings) result(output)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this
  type(PrintSettings),     intent(in), optional :: settings
  type(String), allocatable                     :: output(:)
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  output = this%write()
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end function

recursive function str_StringsWriteables_String(this,separating_line, &
   & settings) result(output)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  type(String),            intent(in), optional :: separating_line
  type(PrintSettings),     intent(in), optional :: settings
  type(String), allocatable                     :: output(:)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i,ialloc
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  allocate(sections(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    sections(i) = StringArray(str(this(i)))
  enddo
  
  output = str(join(sections,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end function

recursive function str_StringsWriteables_character(this,separating_line, &
   & settings) result(output)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  character(*),            intent(in)           :: separating_line
  type(PrintSettings),     intent(in), optional :: settings
  type(String), allocatable                     :: output(:)
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  output = str(this,str(separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end function

! ----------------------------------------------------------------------
! Provides print_lines() for types which extend StringsWriteable.
! ----------------------------------------------------------------------
subroutine print_lines_StringsWriteable(this,settings)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this
  type(PrintSettings),     intent(in), optional :: settings
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end subroutine

subroutine print_lines_StringsWriteables_String(this,separating_line,settings)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  type(String),            intent(in), optional :: separating_line
  type(PrintSettings),     intent(in), optional :: settings
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end subroutine

subroutine print_lines_StringsWriteables_character(this,separating_line, &
   & settings)
  implicit none
  
  class(StringsWriteable), intent(in)           :: this(:)
  character(*),            intent(in)           :: separating_line
  type(PrintSettings),     intent(in), optional :: settings
  
  if (present(settings)) then
    call set_print_settings(settings)
  endif
  
  call print_lines(str(this,separating_line))
  
  if (present(settings)) then
    call unset_print_settings()
  endif
end subroutine
end module

! ======================================================================
! An example module showing how to extend the StringsWriteable type.
! ======================================================================
module caesar_StringsWriteable_example_module
  use caesar_io_basic_module
  
  use caesar_strings_writeable_module
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
