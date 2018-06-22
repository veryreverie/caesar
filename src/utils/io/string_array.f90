! ======================================================================
! An array of type String.
! Exists to allow heterogeneous 2-D arrays of type String.
! ======================================================================
module string_array_submodule
  use io_basic_module
  implicit none
  
  private
  
  public :: StringArray
  public :: size
  public :: str
  public :: split_into_sections
  public :: operator(//)
  public :: join
  
  type :: StringArray
    type(String), allocatable :: strings(:)
  end type
  
  interface StringArray
    module procedure new_StringArray_Strings
  end interface
  
  interface size
    module procedure size_StringArray
  end interface
  
  interface str
    module procedure str_StringArray
  end interface
  
  interface split_into_sections
    module procedure split_into_sections_Strings_character
    module procedure split_into_sections_Strings_String
    module procedure split_into_sections_StringArray_character
    module procedure split_into_sections_StringArray_String
  end interface
  
  interface operator(//)
    module procedure concatenate_StringArray_StringArray
    module procedure concatenate_StringArray_String
    module procedure concatenate_StringArray_character
    module procedure concatenate_StringArray_Strings
    module procedure concatenate_String_StringArray
    module procedure concatenate_character_StringArray
    module procedure concatenate_Strings_StringArray
  end interface
  
  interface join
    module procedure join_StringArrays_String
    module procedure join_StringArrays_character
  end interface
contains

! --------------------------------------------------
! Basic functionality:
!    - Constructor.
!    - The size() function.
!    - The str() function.
! --------------------------------------------------
function new_StringArray_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(StringArray)        :: this
  
  this%strings = input
end function

function size_StringArray(this) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  integer                       :: output
  
  output = size(this%strings)
end function

function str_StringArray(this) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  output = this%strings
end function

! ----------------------------------------------------------------------
! Splits the array into sections.
! Splits by one or more strings which match delimiter.
! Delimiter defaults to an empty string, ''.
! ----------------------------------------------------------------------
function split_into_sections_Strings_character(this,delimiter) result(output)
  implicit none
  
  type(String), intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(StringArray), allocatable     :: output(:)
  
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
    if (this(i)==delimiter_string) then
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
    output(i) = StringArray(this(first_lines(i):last_lines(i)))
  enddo
end function

function split_into_sections_Strings_String(this,delimiter) result(output)
  implicit none
  
  type(String), intent(in)       :: this(:)
  type(String), intent(in)       :: delimiter
  type(StringArray), allocatable :: output(:)
  
  output = split_into_sections(this,char(delimiter))
end function

function split_into_sections_StringArray_character(this,delimiter) &
   & result(output)
  implicit none
  
  type(StringArray), intent(in)           :: this
  character(*),      intent(in), optional :: delimiter
  type(StringArray), allocatable          :: output(:)
  
  output = split_into_sections(this%strings,delimiter)
end function

function split_into_sections_StringArray_String(this,delimiter) result(output)
  implicit none
  
  type(StringArray), intent(in)  :: this
  type(String),      intent(in)  :: delimiter
  type(StringArray), allocatable :: output(:)
  
  output = split_into_sections(this%strings,char(delimiter))
end function

! ----------------------------------------------------------------------
! Concatenate StringArrays and another string-like structure.
! ----------------------------------------------------------------------
function concatenate_StringArray_StringArray(this,that) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  type(StringArray), intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([this%strings, that%strings])
end function

function concatenate_StringArray_String(this,that) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  type(String),      intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([this%strings, that])
end function

function concatenate_StringArray_character(this,that) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  character(*),      intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([this%strings, str(that)])
end function

function concatenate_StringArray_Strings(this,that) result(output)
  implicit none
  
  type(StringArray), intent(in) :: this
  type(String),      intent(in) :: that(:)
  type(StringArray)             :: output
  
  output = StringArray([this%strings, that])
end function

function concatenate_String_StringArray(this,that) result(output)
  implicit none
  
  type(String),      intent(in) :: this
  type(StringArray), intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([this, that%strings])
end function

function concatenate_character_StringArray(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  type(StringArray), intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([str(this), that%strings])
end function

function concatenate_Strings_StringArray(this,that) result(output)
  implicit none
  
  type(String),      intent(in) :: this(:)
  type(StringArray), intent(in) :: that
  type(StringArray)             :: output
  
  output = StringArray([this, that%strings])
end function

function join_StringArrays_String(input,separating_line) result(output)
  implicit none
  
  type(StringArray), intent(in)           :: input(:)
  type(String),      intent(in), optional :: separating_line
  type(StringArray)                       :: output
  
  integer                   :: no_strings
  type(String), allocatable :: strings(:)
  
  integer :: i,j,ialloc
  
  if (size(input)==0) then
    strings = [String::]
  else
    no_strings = size(input(1))
    do i=2,size(input)
      no_strings = no_strings+size(input(i))
      if (present(separating_line)) then
        no_strings = no_strings + 1
      endif
    enddo
    
    allocate(strings(no_strings), stat=ialloc); call err(ialloc)
    
    strings(:size(input(1))) = input(1)%strings
    j = size(input(1))
    do i=2,size(input)
      if (present(separating_line)) then
        strings(j+1) = separating_line
        j = j+1
      endif
      strings(j+1:j+size(input(i))) = input(i)%strings
      j = j+size(input(i))
    enddo
  endif
  
  output = StringArray(strings)
end function

function join_StringArrays_character(input,separating_line) result(output)
  implicit none
  
  type(StringArray), intent(in) :: input(:)
  character(*),      intent(in) :: separating_line
  type(StringArray)             :: output
  
  output = join(input,str(separating_line))
end function
end module
