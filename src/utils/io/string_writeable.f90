! ======================================================================
! An abstract type, which allows extended types to be written to a String.
! ======================================================================
! Any type which extends StringWriteable can be:
!    - converted to String, using str(this).
!    - concatenated, using string//this or character//this.
!    - printed to stdout, using print_line(this).
!    - printed to file, using file%print_line(this).
! An array of StringWriteables can be printed across multiple lines using
!    print_lines(this(:)) or file%print_lines(this(:)),
!    with the optional keyword separating_line to separate elements.
! Any type which extends StringWriteable must overload %write().
! See example module below for how to extend this type.
module string_writeable_submodule
  use io_basic_module
  use abstract_module
  
  use string_array_submodule
  implicit none
  
  private
  
  public :: StringWriteable
  public :: str
  public :: operator(//)
  public :: join
  public :: print_line
  public :: print_lines
  
  type, abstract, extends(NoDefaultConstructor) :: StringWriteable
  contains
    procedure(write_StringWriteable), deferred :: write
  end type
  
  abstract interface
    recursive function write_StringWriteable(this) result(output)
      import String
      import StringWriteable
      implicit none
      
      class(StringWriteable), intent(in) :: this
      type(String)                       :: output
    end function
  end interface
  
  interface str
    module procedure str_StringWriteable
    module procedure str_StringWriteables_character
    module procedure str_StringWriteables_String
  end interface
    
  interface operator(//)
    module procedure concatenate_StringWriteable_character
    module procedure concatenate_character_StringWriteable
    module procedure concatenate_StringWriteable_String
    module procedure concatenate_String_StringWriteable
  end interface
  
  interface join
    module procedure join_StringWriteable
  end interface
  
  interface print_line
    module procedure print_line_StringWriteable
  end interface
  
  interface print_lines
    module procedure print_lines_StringWriteables_character
    module procedure print_lines_StringWriteables_String
  end interface
contains

! ----------------------------------------------------------------------
! The str() function, which converts to string types.
! ----------------------------------------------------------------------
recursive function str_StringWriteable(this) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this
  type(String)                       :: output
  
  output = this%write()
end function

recursive function str_StringWriteables_character(this,separating_line) &
   & result(output)
  implicit none
  
  class(StringWriteable), intent(in)           :: this(:)
  character(*),           intent(in), optional :: separating_line
  type(String), allocatable                    :: output(:)
  
  integer :: i,ialloc
  
  if (present(separating_line)) then
    allocate(output(2*size(this)-1), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(2*i-1) = str(this(i))
      if (i<size(this)) then
        output(2*i) = separating_line
      endif
    enddo
  else
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = str(this(i))
    enddo
  endif
end function

recursive function str_StringWriteables_String(this,separating_line) &
   & result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this(:)
  type(String),           intent(in) :: separating_line
  type(String), allocatable          :: output(:)
  
  output = str(this,char(separating_line))
end function

! ----------------------------------------------------------------------
! Concatenation of string types and StringWriteable types.
! ----------------------------------------------------------------------
! String = StringWriteable//character
recursive function concatenate_StringWriteable_character(this,that) &
   & result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this
  character(*),           intent(in) :: that
  type(String)                       :: output
  
  output = str(this)//that
end function

! String = character//StringWriteable
recursive function concatenate_character_StringWriteable(this,that) &
   & result(output)
  implicit none
  
  character(*),           intent(in) :: this
  class(StringWriteable), intent(in) :: that
  type(String)                       :: output
  
  output = this//str(that)
end function

! String = StringWriteable//String
recursive function concatenate_StringWriteable_String(this,that) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this
  type(String),           intent(in) :: that
  type(String)                       :: output
  
  output = str(this)//that
end function

! String = String//StringWriteable
recursive function concatenate_String_StringWriteable(this,that) result(output)
  implicit none
  
  type(String),           intent(in) :: this
  class(StringWriteable), intent(in) :: that
  type(String)                       :: output
  
  output = this//str(that)
end function

! ----------------------------------------------------------------------
! Convert a StringWriteable array to an array of Strings,
!    then concatenate them into a single string.
! ----------------------------------------------------------------------
recursive function join_StringWriteable(this,delimiter) result(output)
  implicit none
  
  class(StringWriteable), intent(in)           :: this(:)
  character(*),           intent(in), optional :: delimiter
  type(String)                                 :: output
  
  output = join(str(this), delimiter)
end function

! ----------------------------------------------------------------------
! Provides print_line and print_lines for types which extend StringWriteable.
! ----------------------------------------------------------------------
subroutine print_line_StringWriteable(this,indent)
  implicit none
  
  class(StringWriteable), intent(in)           :: this
  integer,                intent(in), optional :: indent
  
  call print_line(str(this),indent)
end subroutine

subroutine print_lines_StringWriteables_character(this,indent,separating_line)
  implicit none
  
  class(StringWriteable), intent(in)           :: this(:)
  integer,                intent(in), optional :: indent
  character(*),           intent(in), optional :: separating_line
  
  call print_lines(str(this,separating_line),indent)
end subroutine

subroutine print_lines_StringWriteables_String(this,indent,separating_line)
  implicit none
  
  class(StringWriteable), intent(in)           :: this(:)
  integer,                intent(in), optional :: indent
  type(String),           intent(in)           :: separating_line
  
  call print_lines(str(this,separating_line),indent)
end subroutine
end module

! ======================================================================
! An example module showing how to extend StringWriteable.
! ======================================================================
module string_writeable_example_submodule
  use io_basic_module
  
  use string_writeable_submodule
  implicit none
  
  private
  
  public :: StringWriteableExample
  
  type, extends(StringWriteable) :: StringWriteableExample
    integer :: contents
  contains
    procedure, public :: write => write_StringWriteableExample
  end type
contains

recursive function write_StringWriteableExample(this) result(output)
  implicit none
  
  class(StringWriteableExample), intent(in) :: this
  type(String)                              :: output
  
  ! Select type needed to call non-polymorphic procedures, and to ensure that
  !    write() is overloaded by any type which extends StringWriteableExample.
  select type(this); type is(StringWriteableExample)
    output = str(this%contents)
  class default
    call print_line(CODE_ERROR//': Called the StringWriteableExample version &
       &of write() from a type other than StringWriteableExample.')
    call err()
  end select
end function
end module
