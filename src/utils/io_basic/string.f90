! ======================================================================
! A simple heap-allocated String class.
! Allows for inhomogeneous character arrays for e.g. storing files.
! ======================================================================
module string_module
  use precision_module
  use string_base_module
  use error_module
  implicit none
  
  private
  
  public :: String
  public :: str
  public :: operator(==)
  public :: operator(/=)
  public :: char
  public :: assignment(=)
  public :: len
  public :: join
  public :: operator(//)
  public :: lower_case
  public :: spaces
  public :: split_line
  public :: slice
  public :: trim
  public :: replace
  
  ! The class itself.
  type, extends(StringBase) :: String
  end type
  
  ! Interfaces
  interface assignment(=)
    module procedure assign_character_String
  end interface
  
  interface str
    module procedure str_character
    module procedure str_String
  end interface
  
  interface len
    module procedure len_String
  end interface
  
  interface join
    module procedure join_String
  end interface
  
  interface operator(//)
    module procedure concatenate_character_String
    module procedure concatenate_String_character
    module procedure concatenate_String_String
  end interface
  
  interface lower_case
    module procedure lower_case_character
    module procedure lower_case_String
  end interface
  
  interface split_line
    module procedure split_line_character
    module procedure split_line_String
  end interface
  
  interface slice
    module procedure slice_character
    module procedure slice_String
  end interface
  
  interface trim
    module procedure trim_String
  end interface
  
  interface replace
    module procedure replace_character
    module procedure replace_String
  end interface
contains

! --------------------------------------------------
! Assignment.
! --------------------------------------------------
! character = String
subroutine assign_character_String(output,input)
  implicit none
  
  character(*), intent(out) :: output
  type(String), intent(in)  :: input
  
  output = char(input)
end subroutine

! ----------------------------------------------------------------------
! Conversion to String
! ----------------------------------------------------------------------
impure elemental function str_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  type(String)             :: output
  
  output = this
end function

impure elemental function str_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Unary operators
! ----------------------------------------------------------------------

! String length. Equivalent to the character len() function.
function len_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = len(char(this))
end function

! Joins a String(:) array into one String.
function join_String(this,delimiter) result(output)
  implicit none
  
  type(String), intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  ! Temporary variables.
  type(String) :: delimiter_character
  integer      :: i
  
  if (present(delimiter)) then
    delimiter_character = delimiter
  else
    delimiter_character = ' '
  endif
  
  if (size(this)==0) then
    output = ''
  else
    output = this(1)
    do i=2,size(this)
      output = output//delimiter_character//this(i)
    enddo
  endif
end function

! ----------------------------------------------------------------------
! Concatenation
! ----------------------------------------------------------------------
                                     
! --------------------------------------------------
! Concatenation of string types and string types.
! --------------------------------------------------

! String = character//String
function concatenate_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = this//char(that)
end function

! String = String//character
function concatenate_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  type(String)              :: output
  
  output = char(this)//that
end function

! String = String//String
function concatenate_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = char(this)//char(that)
end function

! --------------------------------------------------
! Converts a string to lower case
! --------------------------------------------------
impure elemental function lower_case_character(input) result(output)
  implicit none
  
  character(*), intent(in) :: input
  character(len(input))    :: output
  
  character(*), parameter :: lower_chars = "abcdefghijklmnopqrstuvwxyz"
  character(*), parameter :: upper_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
  
  integer :: i,j
  
  output = input
  
  do i=1,len(output)
    j = index(upper_chars, output(i:i))
    if (j/=0) then
      output(i:i) = lower_chars(j:j)
    endif
  enddo
end function

impure elemental function lower_case_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = str(lower_case(char(this)))
end function

! --------------------------------------------------
! Returns a string consisting of the requested number of spaces.
! --------------------------------------------------
function spaces(no_spaces) result(output)
  implicit none
  
  integer, intent(in)       :: no_spaces
  character(:), allocatable :: output
  
  integer :: i,ialloc
  
  allocate(character(no_spaces) :: output, stat=ialloc); call err(ialloc)
  
  do i=1,no_spaces
    output(i:i) = ' '
  enddo
end function

! --------------------------------------------------
! Split a string by a given delimiter.
! --------------------------------------------------
! If no delimiter is specified, splits by whitespace.
function split_line_character(input,delimiter) result(output)
  implicit none
  
  character(*), intent(in)           :: input
  character(1), intent(in), optional :: delimiter
  type(String), allocatable          :: output(:)
  
  ! The position of a delimiter.
  integer :: first
  ! The position of the next delimiter after 'first'.
  integer :: second
  
  ! The number of characters before the next space or tab.
  integer :: next_space
  integer :: next_tab
  
  
  first = 0
  second = 0
  output = [String::]
  do
    ! Search after previously found delimiter.
    first = second
    ! Exit if entire word split.
    if (first == len(input)+1) then
      exit
    endif
    ! Find the next delimiter.
    if (present(delimiter)) then
      second = first + index(input(first+1:), delimiter)
    else
      ! If delimiter is not set, find the next space or tab character.
      next_space = first + index(input(first+1:), ' ')
      next_tab   = first + index(input(first+1:), '	') ! N.B. '[TAB]'
      if (next_space==first) then
        ! No space found.
        second = next_tab
      elseif (next_tab==first) then
        ! No tab found.
        second = next_space
      else
        second = min(next_tab,next_space)
      endif
    endif
    ! If second==first there is no next delimiter. Parse the final token.
    if (second == first) then
      second = len(input)+1
    endif
    ! If second==first+1, there are multiple delimiters in a row.
    ! They are treated as a single delimiter.
    if (second == first+1) then
      cycle
    endif
    ! Append the token to the output.
    output = [output, str(input(first+1:second-1))]
  enddo
end function

function split_line_String(input,delimiter) result(output)
  implicit none
  
  type(String), intent(in)           :: input
  character(1), intent(in), optional :: delimiter
  type(String), allocatable          :: output(:)
  
  output = split_line(char(input),delimiter)
end function

! --------------------------------------------------
! Takes a slice of a String. slice(String,first,last) = character(first:last).
! --------------------------------------------------
function slice_character(input,first,last) result(output)
  implicit none
  
  character(*), intent(in) :: input
  integer,      intent(in) :: first
  integer,      intent(in) :: last
  type(String)             :: output
  
  output = input(first:last)
end function

impure elemental function slice_String(input,first,last) result(output)
  implicit none
  
  type(String), intent(in) :: input
  integer,      intent(in) :: first
  integer,      intent(in) :: last
  type(String)             :: output
  
  output = slice(char(input),first,last)
end function

! --------------------------------------------------
! Removes trailing spaces.
! --------------------------------------------------
impure elemental function trim_String(input) result(output)
  implicit none
  
  type(String), intent(in) :: input
  type(String)             :: output
  
  output = trim(adjustl(char(input)))
end function

! --------------------------------------------------
! Finds and replaces one character with another.
! --------------------------------------------------
function replace_character(input,from,to) result(output)
  implicit none
  
  character(*), intent(in) :: input
  character(1), intent(in) :: from
  character(1), intent(in) :: to
  type(String)             :: output
  
  integer :: i
  
  output = ''
  do i=1,len(input)
    if (input(i:i)==from) then
      output = output//to
    else
      output = output//input(i:i)
    endif
  enddo
end function

function replace_String(input,from,to) result(output)
  implicit none
  
  type(String), intent(in) :: input
  character(*), intent(in) :: from
  character(*), intent(in) :: to
  type(String)             :: output
  
  output = replace(char(input),from,to)
end function
end module
