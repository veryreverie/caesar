! ======================================================================
! An array of type String.
! Exists to allow heterogeneous 2-D arrays of type String.
! ======================================================================
module string_array_submodule
  use error_submodule
  use string_submodule
  use printable_submodule
  implicit none
  
  private
  
  public :: StringArray
  public :: size
  public :: split
  
  type, extends(Printable) :: StringArray
    type(String), allocatable :: strings(:)
  contains
    procedure, public :: to_String => to_String_StringArray
  end type
  
  interface StringArray
    module procedure new_StringArray_Strings
  end interface
  
  interface size
    module procedure size_StringArray
  end interface
  
  interface split
    module procedure split_StringArray_character
    module procedure split_StringArray_String
  end interface
contains

! --------------------------------------------------
! Basic functionality:
!    - Converting between String(:) and StringArray.
!    - Assignments between String(:) and StringArray.
!    - The size() function.
! --------------------------------------------------
function new_StringArray_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StringArray)        :: this
  
  this%strings = input
end function

function to_String_StringArray(this) result(output)
  implicit none
  
  class(StringArray), intent(in) :: this
  type(String), allocatable      :: output(:)
  
  output = this%strings
end function

function size_StringArray(this) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  integer                       :: output
  
  output = size(this%strings)
end function

! ----------------------------------------------------------------------
! Splits the array into sections.
! Splits by one or more strings which match delimiter.
! Delimiter defaults to an empty string, ''.
! ----------------------------------------------------------------------
function split_StringArray_character(this,delimiter) result(output)
  implicit none
  
  type(StringArray), intent(in)           :: this
  character(*),      intent(in), optional :: delimiter
  type(StringArray), allocatable          :: output(:)
  
  type(String) :: delimiter_string
  
  integer :: no_sections
  logical :: reading_section
  
  integer, allocatable :: first_lines(:)
  integer, allocatable :: last_lines(:)
  
  integer :: i,ialloc
  
  if (present(delimiter)) then
    delimiter_string = delimiter
  else
    delimiter_string = ''
  endif
  
  allocate( first_lines(size(this)), &
          & last_lines(size(this)),  &
          & stat=ialloc); call err(ialloc)
  no_sections = 0
  reading_section = .false.
  do i=1,size(this)
    if (this%strings(i)==delimiter_string) then
      ! This line is the delimiter string.
      ! If reading a section, then the end of that section is the line above.
      if (reading_section) then
        last_lines(no_sections) = i-1
        reading_section = .false.
      endif
    else
      ! This line is not the delimiter string.
      ! If not reading a section, then this line is the start of a new section.
      if (.not. reading_section) then
        no_sections = no_sections+1
        first_lines(no_sections) = i
        reading_section = .true.
      endif
    endif
  enddo
  
  ! If a section is still being read, then that section ends on the last line
  !    of the file.
  if (reading_section) then
    last_lines(no_sections) = size(this)
  endif
  
  allocate(output(no_sections), stat=ialloc); call err(ialloc)
  do i=1,no_sections
    output(i) = StringArray(this%strings(first_lines(i):last_lines(i)))
  enddo
end function

function split_StringArray_String(this,delimiter) result(output)
  implicit none
  
  type(StringArray), intent(in)  :: this
  type(String),      intent(in)  :: delimiter
  type(StringArray), allocatable :: output(:)
  
  output = split(this,char(delimiter))
end function
end module
