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
  use string_submodule
  use string_array_submodule
  use error_submodule
  use print_submodule
  implicit none
  
  private
  
  public :: StringWriteable
  !public :: assignment(=)
  public :: operator(//)
  public :: str
  public :: print_line
  public :: print_lines
  
  type, abstract :: StringWriteable
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
  
  !interface assignment(=)
  !  module procedure assign_String_StringWriteable
  !  module procedure assign_Strings_StringWriteables
  !  module procedure assign_StringArray_StringWriteables
  !end interface
    
  interface operator(//)
    module procedure concatenate_StringWriteable_character
    module procedure concatenate_character_StringWriteable
    module procedure concatenate_StringWriteable_String
    module procedure concatenate_String_StringWriteable
  end interface
  
  interface str
    module procedure str_StringWriteable_0d
    module procedure str_StringWriteable_1d
    module procedure str_StringWriteable_2d
  end interface
  
  interface print_line
    module procedure print_line_StringWriteable
  end interface
  
  interface print_lines
    module procedure print_lines_StringWriteables_character
    module procedure print_lines_StringWriteables_String
  end interface
contains

!! ----------------------------------------------------------------------
!! Assign a string type from a StringWriteable type.
!! ----------------------------------------------------------------------
!recursive subroutine assign_String_StringWriteable(output,input)
!  implicit none
!  
!  type(String),           intent(out) :: output
!  class(StringWriteable), intent(in)  :: input
!  
!  output = input%write()
!end subroutine
!
!recursive subroutine assign_Strings_StringWriteables(output,input)
!  implicit none
!  
!  type(String), allocatable, intent(out) :: output(:)
!  class(StringWriteable),    intent(in)  :: input(:)
!  
!  integer :: i,ialloc
!  
!  allocate(output(size(input)), stat=ialloc); call err(ialloc)
!  do i=1,size(input)
!    output(i) = input(i)
!  enddo
!end subroutine
!
!recursive subroutine assign_StringArray_StringWriteables(output,input)
!  implicit none
!  
!  type(StringArray),      intent(out) :: output
!  class(StringWriteable), intent(in)  :: input(:)
!  
!  output%strings = input
!end subroutine

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
  
  output = this%write()//that
end function

! String = character//StringWriteable
recursive function concatenate_character_StringWriteable(this,that) &
   & result(output)
  implicit none
  
  character(*),           intent(in) :: this
  class(StringWriteable), intent(in) :: that
  type(String)                       :: output
  
  output = this//that%write()
end function

! String = StringWriteable//String
recursive function concatenate_StringWriteable_String(this,that) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this
  type(String),           intent(in) :: that
  type(String)                       :: output
  
  output = this%write()//that
end function

! String = String//StringWriteable
recursive function concatenate_String_StringWriteable(this,that) result(output)
  implicit none
  
  type(String),           intent(in) :: this
  class(StringWriteable), intent(in) :: that
  type(String)                       :: output
  
  output = this//that%write()
end function

! ----------------------------------------------------------------------
! The str() function, which converts to string types.
! ----------------------------------------------------------------------
! N.B. can't use impure elemental because this must be recursive.
recursive function str_StringWriteable_0d(this) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this
  type(String)                       :: output
  
  output = this%write()
end function

recursive function str_StringWriteable_1d(this) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this(:)
  type(String), allocatable          :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = str(this(i))
  enddo
end function

recursive function str_StringWriteable_2d(this) result(output)
  implicit none
  
  class(StringWriteable), intent(in) :: this(:,:)
  type(String), allocatable          :: output(:,:)
  
  integer :: i,ialloc
  
  allocate(output(size(this,1),size(this,2)), stat=ialloc); call err(ialloc)
  do i=1,size(this,2)
    output(:,i) = str(this(:,i))
  enddo
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
  
  integer :: i
  
  do i=1,size(this)
    call print_line(this(i), indent)
    if (present(separating_line)) then
      call print_line(separating_line, indent)
    endif
  enddo
end subroutine

subroutine print_lines_StringWriteables_String(this,indent,separating_line)
  implicit none
  
  class(StringWriteable), intent(in)           :: this(:)
  integer,                intent(in), optional :: indent
  type(String),           intent(in)           :: separating_line
  
  call print_lines(this, indent, char(separating_line))
end subroutine
end module

! ======================================================================
! An example module showing how to extend StringWriteable.
! ======================================================================
module string_writeable_example_submodule
  use string_submodule
  use string_writeable_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: StringWriteableExample
  
  type, extends(StringWriteable) :: StringWriteableExample
    character(1) :: contents
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
