! ======================================================================
! A simple heap-allocated String class
! ======================================================================
module string_module
  use constants_module, only : dp
  implicit none
  
  private
  
  ! ----------------------------------------------------------------------
  ! Public interface
  ! ----------------------------------------------------------------------
  
  public :: Stringable
  public :: String
  
  ! Conversions between classes
  public :: str   ! Conversion to String.
  public :: char  ! Conversion from String to character.
  
  public :: assignment(=)
  
  ! Unary operators.
  public :: len            ! Character-like len(String).
  public :: lower_case     ! Convert to lower case.
  public :: split          ! Split into String(:), by default by spaces.
  public :: join           ! Join into single String, by default with spaces.
  public :: pad_int_to_str ! left pads integers without '-' signs with a ' '
  public :: trim           ! Removes trailing spaces.
  public :: slice          ! slice(String,first,last)=character(first:last).
  
  ! Concatenate to String.
  public :: operator(//)
  
  ! ----------------------------------------------------------------------
  ! The String class.
  ! ----------------------------------------------------------------------
  
  ! The String class, containing an allocatable character string.
  ! Allows for inhomogeneous character arrays for e.g. storing files.
  type :: String
    character(:), allocatable, private :: contents_
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
  
  ! ----------------------------------------------------------------------
  ! The Stringable class.
  ! ----------------------------------------------------------------------
  
  ! An abstract type, which allows extended types to be turned into strings.
  ! Any type which extends Stringable can be:
  !    - converted to String, using string=this or str(this).
  !    - concatenated, using string//this or character//this.
  !    - printed to stdout, using print_line(this).
  !    - printed to file, using file%print_line(this).
  ! See example module below for how to use this module.
  type, abstract :: Stringable
  contains
    procedure(str_Stringable_2), deferred :: str
  end type
  
  abstract interface
    function str_Stringable_2(this) result(output)
      import String
      import Stringable
      implicit none
      
      class(Stringable), intent(in) :: this
      type(String)                  :: output
    end function
  end interface
  
  ! ----------------------------------------------------------------------
  ! Interfaces
  ! ----------------------------------------------------------------------
  interface assignment(=)
    module procedure assign_String_character
    module procedure assign_String_String
    module procedure assign_String_Stringable
    module procedure assign_String_logical
    module procedure assign_String_integer
    module procedure assign_String_real
    module procedure assign_String_complex
  end interface
  
  interface assignment(=)
    module procedure assign_character_String
  end interface
  
  interface str
    module procedure str_character
    module procedure str_String
    module procedure str_Stringable
    module procedure str_integer
    module procedure str_real
    module procedure str_logical
    module procedure str_complex
  end interface
  
  interface char
    module procedure char_String
  end interface
  
  interface len
    module procedure len_String
  end interface
  
  interface lower_case
    module procedure lower_case_character
    module procedure lower_case_String
  end interface
  
  interface split
    module procedure split_character
    module procedure split_String
  end interface
  
  interface join
    module procedure join_String
    module procedure join_real
    module procedure join_integer
    module procedure join_logical
    module procedure join_complex
  end interface
  
  interface trim
    module procedure trim_String
  end interface
  
  interface slice
    module procedure slice_character
    module procedure slice_String
  end interface
  
  interface operator(//)
    module procedure concatenate_character_String
    module procedure concatenate_String_character
    module procedure concatenate_String_String
    
    module procedure concatenate_character_Stringable
    module procedure concatenate_Stringable_character
    module procedure concatenate_String_Stringable
    module procedure concatenate_Stringable_String
                                     
    module procedure concatenate_character_integer
    module procedure concatenate_integer_character
    module procedure concatenate_String_integer
    module procedure concatenate_integer_String
    
    module procedure concatenate_character_real
    module procedure concatenate_real_character
    module procedure concatenate_String_real
    module procedure concatenate_real_String
    
    module procedure concatenate_character_logical
    module procedure concatenate_logical_character
    module procedure concatenate_String_logical
    module procedure concatenate_logical_String
    
    module procedure concatenate_String_complex
    module procedure concatenate_complex_String
    module procedure concatenate_character_complex
    module procedure concatenate_complex_character
    
    module procedure concatenate_character_integers
    module procedure concatenate_integers_character
    module procedure concatenate_String_integers
    module procedure concatenate_integers_String
    
    module procedure concatenate_character_reals
    module procedure concatenate_reals_character
    module procedure concatenate_String_reals
    module procedure concatenate_reals_String
    
    module procedure concatenate_character_logicals
    module procedure concatenate_logicals_character
    module procedure concatenate_String_logicals
    module procedure concatenate_logicals_String
    
    module procedure concatenate_character_complexes
    module procedure concatenate_complexes_character
    module procedure concatenate_String_complexes
    module procedure concatenate_complexes_String
  end interface
contains

! ----------------------------------------------------------------------
! String operations involving the private variable contents_
! ----------------------------------------------------------------------
! Assignment.
! String = character
pure subroutine assign_String_character(output,input)
  implicit none
  
  type(String), intent(out) :: output
  character(*), intent(in)  :: input
  
  output%contents_ = input
end subroutine

! String length. Equivalent to the character len() function.
pure function len_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  if (allocated(this%contents_)) then
    output = len(this%contents_)
  else
    output = 0
  endif
end function

! Converts a String to a character(*). Effectively a getter for contents_.
pure function char_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  character(len(this))     :: output
  
  if (allocated(this%contents_)) then
    output = this%contents_
  else
    output = ''
  endif
end function

! ----------------------------------------------------------------------
! String operations not involving the private variable contents_
! ----------------------------------------------------------------------
! N.B. the number of procedures accessing contents_ directly is intentionally
!    limited for stability reasons.
! The above procedures behave well if the String has not been allocated,
!    and this good behaviour is automatically passed to the procedures below.

! --------------------------------------------------
! Assignment.
! --------------------------------------------------
! String = String
pure subroutine assign_String_String(output,input)
  implicit none
  
  type(String),  intent(out) :: output
  class(String), intent(in)  :: input
  
  output = char(input)
end subroutine

! String = Stringable
subroutine assign_String_Stringable(output,input)
  implicit none
  
  type(String),      intent(out) :: output
  class(Stringable), intent(in)  :: input
  
  output = input%str()
end subroutine

! String = logical
pure subroutine assign_String_logical(output,input)
  implicit none
  
  type(String), intent(out) :: output
  logical,      intent(in)  :: input
  
  if (input) then
    output = 'T'
  else
    output = 'F'
  endif
end subroutine

! String = integer
pure subroutine assign_String_integer(output,input)
  implicit none
  
  type(String), intent(out) :: output
  integer,      intent(in)  :: input
  
  integer, parameter   :: int_width = 12
  character(int_width) :: int_string
  
  write(int_string,"(I0)") input
  output = trim(int_string)
end subroutine

! String = real
pure subroutine assign_String_real(output,input)
  implicit none
  
  type(String), intent(out) :: output
  real(dp),     intent(in)  :: input
  
  integer, parameter    :: real_width = 25
  integer, parameter    :: decimal_places = 17
  type(String)          :: format_string
  character(real_width) :: real_string
  
  format_string = "(ES"//real_width//'.'//decimal_places//")"
  write(real_string, char(format_string)) input
  
  output = real_string
end subroutine

! String = complex
pure subroutine assign_String_complex(output,input)
  implicit none
  
  type(String), intent(out) :: output
  complex(dp),  intent(in)  :: input
  
  type(String) :: real_string
  type(String) :: imag_string
  
  real_string = real(input)
  imag_string = aimag(input)
  
  imag_string = trim(imag_string)
  
  if (slice(imag_string,1,1)/='-') then
    imag_string = '+'//imag_string
  endif
  
  output = real_string//imag_string//'i'
end subroutine

! character = String
pure subroutine assign_character_String(output,input)
  implicit none
  
  character(*), intent(out) :: output
  type(String), intent(in)  :: input
  
  output = char(input)
end subroutine

! ----------------------------------------------------------------------
! Conversion to String
! ----------------------------------------------------------------------
elemental function str_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  type(String)             :: output
  
  output = this
end function

elemental function str_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = this
end function

impure elemental function str_Stringable(this) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String)                  :: output
  
  output = this
end function

elemental function str_logical(this) result(output)
  implicit none
  
  logical, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

elemental function str_integer(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

elemental function str_real(this) result(output)
  implicit none
  
  real(dp), intent(in) :: this
  type(String)         :: output
  
  output = this
end function

elemental function str_complex(this) result(output)
  implicit none
  
  complex(dp), intent(in) :: this
  type(String)            :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Equality
! ----------------------------------------------------------------------
! String==String
elemental function equality_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = char(this)==char(that)
end function

! String==character
elemental function equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = char(this)==that
end function

! character==String
elemental function equality_character_String(this,that) result(output)
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
elemental function non_equality_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
end function

! String/=character
elemental function non_equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
end function

! character/=String
elemental function non_equality_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! Unary operators
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Converts a string to lower case
! ----------------------------------------------------------------------
elemental function lower_case_character(input) result(output)
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

! String = lower_case(String)
elemental function lower_case_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = str(lower_case(char(this)))
end function

pure function split_character(this,delimiter_in) result(output)
  implicit none
  
  character(*),           intent(in) :: this
  character(1), optional, intent(in) :: delimiter_in
  type(String), allocatable          :: output(:)
  
  ! Working variables
  character(1) :: delimiter
  integer      :: first  ! The position of a delimiter
  integer      :: second ! The position of the next delimiter after 'first'
  integer      :: count  ! The number of tokens
  
  if (present(delimiter_in)) then
    delimiter = delimiter_in
  else
    delimiter = ' '
  endif
  
  ! Count the number of tokens in the string
  first = 0
  second = 0
  count = 0
  do
    first = second ! Search after previously found delimiter
    if (first == len(this)+1) exit ! Exit if entire word split
    second = first+index(this(first+1:),delimiter) ! Find the next delimiter
    if (second == first) second = len(this)+1 ! Split the final token
    if (second == first+1) cycle ! Ignore multiple delimiters in a row
    count = count + 1
  enddo
  
  ! Allocate output
  allocate(output(count))
  
  ! Split string
  first = 0
  second = 0
  count = 0
  do
    first = second
    if (first == len(this)+1) exit
    second = first + index(this(first+1:),delimiter)
    if (second == first) second = len(this)+1
    if (second == first+1) cycle
    count = count + 1
    output(count) = this(first+1:second-1)
  enddo
end function

pure function split_String(this,delimiter_in) result(output)
  implicit none
  
  type(String),           intent(in) :: this
  character(1), optional, intent(in) :: delimiter_in
  type(String), allocatable          :: output(:)
  
  if (present(delimiter_in)) then
    output = split(char(this),delimiter_in)
  else
    output = split(char(this))
  endif
end function

! Joins a String(:) array into one String.
pure function join_String(this,delimiter_in) result(output)
  implicit none
  
  type(String), intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter_in
  type(String)                       :: output
  
  ! Temporary variables.
  type(String) :: delimiter
  integer      :: i
  
  if (present(delimiter_in)) then
    delimiter = delimiter_in
  else
    delimiter = ' '
  endif
  
  if (size(this)==0) then
    output = ''
  else
    output = this(1)
    do i=2,size(this)
      output = output//delimiter//this(i)
    enddo
  endif
end function

! Converts real(:) to String(:) and then joins.
pure function join_real(this,delimiter) result(output)
  implicit none
  
  real(dp),     intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  if (present(delimiter)) then
    output = join(str(this),delimiter)
  else
    output = join(str(this))
  endif
end function

! Converts integer(:) to String(:), pads +ve numbers with a space, and joins.
pure function join_integer(this,delimiter) result(output)
  implicit none
  
  integer,      intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  if (present(delimiter)) then
    output = join(pad_int_to_str(this),delimiter)
  else
    output = join(pad_int_to_str(this))
  endif
end function

! Converts logical(:) to String(:) and then joins.
pure function join_logical(this,delimiter) result(output)
  implicit none
  
  logical,      intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  if (present(delimiter)) then
    output = join(str(this),delimiter)
  else
    output = join(str(this))
  endif
end function

! Converts complex(:) to String(:) and then joins.
pure function join_complex(this,delimiter) result(output)
  implicit none
  
  complex(dp),  intent(in)           :: this(:)
  character(*), intent(in), optional :: delimiter
  type(String)                       :: output
  
  if (present(delimiter)) then
    output = join(str(this),delimiter)
  else
    output = join(str(this))
  endif
end function

! Pads an integer with a space to match '-' length.
elemental function pad_int_to_str(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
  if (slice(output,1,1)/='-') then
    output = ' '//output
  endif
end function

! Removes trailing spaces.
elemental function trim_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  type(String)             :: output
  
  output = trim(adjustl(char(this)))
end function

! Takes a slice of a String. slice(String,first,last) = character(first:last).
pure function slice_character(this,first,last) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  integer,       intent(in) :: first
  integer,       intent(in) :: last
  type(String)              :: output
  
  output = this(first:last)
end function

pure function slice_String(this,first,last) result(output)
  implicit none
  
  type(String),  intent(in) :: this
  integer,       intent(in) :: first
  integer,       intent(in) :: last
  type(String)              :: output
  
  output = slice(char(this),first,last)
end function

! ----------------------------------------------------------------------
! Concatenation
! ----------------------------------------------------------------------
                                     
! --------------------------------------------------
! Concatenation of string types and string types.
! --------------------------------------------------

! String = character//String
pure function concatenate_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = this//char(that)
end function

! String = String//character
pure function concatenate_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  type(String)              :: output
  
  output = char(this)//that
end function

! String = String//String
pure function concatenate_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = char(this)//char(that)
end function

! --------------------------------------------------
! Concatenation of string types and Stringable types.
! --------------------------------------------------

! String = character//Stringable
function concatenate_character_Stringable(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//that%str()
end function

! String = Stringable//character
function concatenate_Stringable_character(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  character(*),      intent(in) :: that
  type(String)                  :: output
  
  output = this%str()//that
end function

! String = String//Stringable
function concatenate_String_Stringable(this,that) result(output)
  implicit none
  
  type(String),      intent(in) :: this
  class(Stringable), intent(in) :: that
  type(String)                  :: output
  
  output = this//that%str()
end function

! String = Stringable//String
function concatenate_Stringable_String(this,that) result(output)
  implicit none
  
  class(Stringable), intent(in) :: this
  type(String),      intent(in) :: that
  type(String)                  :: output
  
  output = this%str()//that
end function

! --------------------------------------------------
! Concatenation of string types and default types.
! --------------------------------------------------

! String = String//integer
pure function concatenate_String_integer(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  integer,       intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = integer//String
pure function concatenate_integer_String(this,that) result(output)
  implicit none
  
  integer,       intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//real(dp)
pure function concatenate_String_real(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  real(dp),      intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = real(dp)//String
pure function concatenate_real_String(this,that) result(output)
  implicit none
  
  real(dp),      intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//logical
pure function concatenate_String_logical(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  logical,       intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = logical//String
pure function concatenate_logical_String(this,that) result(output)
  implicit none
  
  logical,       intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//complex
pure function concatenate_String_complex(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  complex(dp),   intent(in) :: that
  type(String)              :: output
  
  output = this//str(that)
end function

! String = complex//String
pure function concatenate_complex_String(this,that) result(output)
  implicit none
  
  complex(dp),   intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = str(this)//that
end function

! String = String//integer(:)
pure function concatenate_String_integers(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  integer,       intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = integer(:)//String
pure function concatenate_integers_String(this,that) result(output)
  implicit none
  
  integer,       intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//real(dp)(:)
pure function concatenate_String_reals(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  real(dp),      intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = real(dp)(:)//String
pure function concatenate_reals_String(this,that) result(output)
  implicit none
  
  real(dp),      intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//logical(:)
pure function concatenate_String_logicals(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  logical,       intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = logical(:)//String
pure function concatenate_logicals_String(this,that) result(output)
  implicit none
  
  logical,       intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = String//complex(:)
pure function concatenate_String_complexes(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  complex(dp),   intent(in) :: that(:)
  type(String)              :: output
  
  output = this//join(that)
end function

! String = complex(:)//String
pure function concatenate_complexes_String(this,that) result(output)
  implicit none
  
  complex(dp),   intent(in) :: this(:)
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = join(this)//that
end function

! String = character//integer
pure function concatenate_character_integer(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer,      intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = integer//character
pure function concatenate_integer_character(this,that) result(output)
  implicit none
  
  integer,      intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//real(dp)
pure function concatenate_character_real(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp),     intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = real(dp)//character
pure function concatenate_real_character(this,that) result(output)
  implicit none
  
  real(dp),     intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//logical
pure function concatenate_character_logical(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical,      intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = logical//character
pure function concatenate_logical_character(this,that) result(output)
  implicit none
  
  logical,      intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//complex
pure function concatenate_character_complex(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp),  intent(in) :: that
  type(String)             :: output
  
  output = this//str(that)
end function

! String = complex//character
pure function concatenate_complex_character(this,that) result(output)
  implicit none
  
  complex(dp),  intent(in) :: this
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = str(this)//that
end function

! String = character//integer(:)
pure function concatenate_character_integers(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  integer,      intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = integer(:)//character
pure function concatenate_integers_character(this,that) result(output)
  implicit none
  
  integer,      intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//real(dp)(:)
pure function concatenate_character_reals(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  real(dp),     intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = real(dp)(:)//character
pure function concatenate_reals_character(this,that) result(output)
  implicit none
  
  real(dp),     intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//logical(:)
pure function concatenate_character_logicals(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  logical,      intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = logical(:)//character
pure function concatenate_logicals_character(this,that) result(output)
  implicit none
  
  logical,      intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function

! String = character//complex(:)
pure function concatenate_character_complexes(this,that) result(output)
  implicit none
  
  character(*), intent(in) :: this
  complex(dp),  intent(in) :: that(:)
  type(String)             :: output
  
  output = this//join(that)
end function

! String = complex(:)//character
pure function concatenate_complexes_character(this,that) result(output)
  implicit none
  
  complex(dp),  intent(in) :: this(:)
  character(*), intent(in) :: that
  type(String)             :: output
  
  output = join(this)//that
end function
end module
