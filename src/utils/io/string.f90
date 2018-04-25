! ======================================================================
! A simple heap-allocated String class.
! Allows for inhomogeneous character arrays for e.g. storing files.
! ======================================================================
module string_submodule
  use precision_module
  use string_base_submodule
  use error_submodule
  implicit none
  
  private
  
  public :: String
  public :: str
  public :: char
  public :: assignment(=)
  public :: len
  public :: join
  public :: operator(//)
  public :: lower_case
  public :: spaces
  public :: split
  public :: slice
  public :: trim
  
  ! The class itself.
  type, extends(StringBase) :: String
  contains
    generic, public :: operator(==) => equality_String_String,    &
                                     & equality_String_character, &
                                     & equality_character_String
    
    generic, public :: operator(/=) => non_equality_String_String,    &
                                     & non_equality_String_character, &
                                     & non_equality_character_String
    
    
    procedure, private             :: equality_String_String
    procedure, private             :: equality_String_character
    procedure, private, pass(that) :: equality_character_String
    
    procedure, private             :: non_equality_String_String
    procedure, private             :: non_equality_String_character
    procedure, private, pass(that) :: non_equality_character_String
  end type
  
  ! Interfaces
  interface assignment(=)
    module procedure assign_String_String
  end interface
  
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
  
  interface split
    module procedure split_character
    module procedure split_String
  end interface
  
  interface slice
    module procedure slice_character
    module procedure slice_String
  end interface
  
  interface trim
    module procedure trim_String
  end interface
contains

! --------------------------------------------------
! Assignment.
! --------------------------------------------------
! String = String
subroutine assign_String_String(output,input)
  implicit none
  
  type(String),  intent(out) :: output
  class(String), intent(in)  :: input
  
  output = char(input)
end subroutine

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
! Equality
! ----------------------------------------------------------------------
! String==String
impure elemental function equality_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = char(this)==char(that)
end function

! String==character
impure elemental function equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = char(this)==that
end function

! character==String
impure elemental function equality_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = this==char(that)
end function

! ----------------------------------------------------------------------
! Non-equality
! ----------------------------------------------------------------------
! String/=String
impure elemental function non_equality_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
end function

! String/=character
impure elemental function non_equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
end function

! character/=String
impure elemental function non_equality_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
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
function split_character(this,delimiter) result(output)
  implicit none
  
  character(*), intent(in)           :: this
  character(1), intent(in), optional :: delimiter
  type(String), allocatable          :: output(:)
  
  ! Set to delimiter if present, or ' ' if not.
  character(1) :: delimiter_character
  
  ! Locations to help parsing.
  integer :: first  ! The position of a delimiter.
  integer :: second ! The position of the next delimiter after 'first'.
  integer :: count  ! The number of tokens.
  
  integer :: ialloc
  
  if (present(delimiter)) then
    delimiter_character = delimiter
  else
    delimiter_character = ' '
  endif
  
  ! Count the number of tokens in the string.
  first = 0
  second = 0
  count = 0
  do
    ! Search after previously found delimiter.
    first = second
    ! Exit if entire word split.
    if (first == len(this)+1) then
      exit
    endif
    ! Find the next delimiter.
    second = first + index(this(first+1:),delimiter_character)
    ! Split the final token.
    if (second == first) then
      second = len(this)+1
    endif
    ! Ignore multiple delimiters in a row.
    if (second == first+1) then
      cycle
    endif
    count = count + 1
  enddo
  
  ! Allocate output.
  allocate(output(count), stat=ialloc); call err(ialloc)
  
  ! Split string. Logic as above, but tokens are transferred to output
  !    rather than just being counted.
  first = 0
  second = 0
  count = 0
  do
    first = second
    if (first == len(this)+1) then
      exit
    endif
    second = first + index(this(first+1:),delimiter_character)
    if (second == first) then
      second = len(this)+1
    endif
    if (second == first+1) then
      cycle
    endif
    count = count + 1
    output(count) = this(first+1:second-1)
  enddo
end function

function split_String(this,delimiter) result(output)
  implicit none
  
  type(String), intent(in)           :: this
  character(1), intent(in), optional :: delimiter
  type(String), allocatable          :: output(:)
  
  output = split(char(this),delimiter)
end function

! --------------------------------------------------
! Takes a slice of a String. slice(String,first,last) = character(first:last).
! --------------------------------------------------
function slice_character(this,first,last) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer,      intent(in) :: first
  integer,      intent(in) :: last
  type(String)             :: output
  
  output = this(first:last)
end function

function slice_String(this,first,last) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer,      intent(in) :: first
  integer,      intent(in) :: last
  type(String)             :: output
  
  output = slice(char(this),first,last)
end function

! --------------------------------------------------
! Removes trailing spaces.
! --------------------------------------------------
impure elemental function trim_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = trim(adjustl(char(this)))
end function
end module
