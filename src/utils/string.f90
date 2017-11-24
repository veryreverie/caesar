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
  
  public :: String ! String class.
  
  ! Conversions between classes
  public :: str   ! Conversion to String.
  public :: char  ! Conversion from String to character.
  public :: int   ! Conversion from String to integer.
  public :: dble  ! Conversion from String to real(dp).
  public :: cmplx ! Conversion from String to complex(dp).
  public :: lgcl  ! Conversion from String to logical.
  
  ! Concatenate to String.
  public :: operator(//)
  
  ! Unary operators.
  public :: len            ! Character-like len(String).
  public :: lower_case     ! Convert to lower case.
  public :: split          ! Split into String(:), by default by spaces.
  public :: join           ! Join into single String, by default with spaces.
  public :: pad_int_to_str ! left pads integers without '-' signs with a ' '
  public :: trim           ! Removes trailing spaces.
  
  ! String slice.
  public :: slice ! slice(String,first,last) = character(first:last).
  
  ! String type.
  type :: String
    character(:), allocatable, private :: contents
  contains
    generic, public :: assignment(=) => assign_String_character, &
                                      & assign_String_String,    &
                                      & assign_String_integer,   &
                                      & assign_String_real,      &
                                      & assign_String_logical,   &
                                      & assign_String_complex,   &
                                      & assign_character_String
    
    generic, public :: operator(//) => concatenate_String_String,    &
                                     & concatenate_String_character, &
                                     & concatenate_character_String, &
                                     & concatenate_String_integer,   &
                                     & concatenate_integer_String,   &
                                     & concatenate_String_real,      &
                                     & concatenate_real_String,      &
                                     & concatenate_String_logical,   &
                                     & concatenate_logical_String,   &
                                     & concatenate_String_complex,   &
                                     & concatenate_complex_String,   &
                                     & concatenate_String_integers,  &
                                     & concatenate_integers_String,  &
                                     & concatenate_String_reals,     &
                                     & concatenate_reals_String,     &
                                     & concatenate_String_logicals,  &
                                     & concatenate_logicals_String,  &
                                     & concatenate_String_complexes, &
                                     & concatenate_complexes_String
    
    generic, public :: operator(==) => equality_String_String,    &
                                     & equality_String_character, &
                                     & equality_character_String
    
    generic, public :: operator(/=) => non_equality_String_String,    &
                                     & non_equality_String_character, &
                                     & non_equality_character_String
    
    procedure, private             :: assign_String_character
    procedure, private             :: assign_String_String
    procedure, private             :: assign_String_integer
    procedure, private             :: assign_String_real
    procedure, private             :: assign_String_logical
    procedure, private             :: assign_String_complex
    procedure, private, pass(that) :: assign_character_String
    
    procedure, private             :: concatenate_String_String
    procedure, private             :: concatenate_String_character
    procedure, private, pass(that) :: concatenate_character_String
    procedure, private             :: concatenate_String_integer
    procedure, private, pass(that) :: concatenate_integer_String
    procedure, private             :: concatenate_String_real
    procedure, private, pass(that) :: concatenate_real_String
    procedure, private             :: concatenate_String_logical
    procedure, private, pass(that) :: concatenate_logical_String
    procedure, private             :: concatenate_String_complex
    procedure, private, pass(that) :: concatenate_complex_String
    procedure, private             :: concatenate_String_integers
    procedure, private, pass(that) :: concatenate_integers_String
    procedure, private             :: concatenate_String_reals
    procedure, private, pass(that) :: concatenate_reals_String
    procedure, private             :: concatenate_String_logicals
    procedure, private, pass(that) :: concatenate_logicals_String
    procedure, private             :: concatenate_String_complexes
    procedure, private, pass(that) :: concatenate_complexes_String
    
    procedure, private             :: equality_String_String
    procedure, private             :: equality_String_character
    procedure, private, pass(that) :: equality_character_String
    
    procedure, private             :: non_equality_String_String
    procedure, private             :: non_equality_String_character
    procedure, private, pass(that) :: non_equality_character_String
  end type
  
  ! ----------------------------------------------------------------------
  ! Interfaces
  ! ----------------------------------------------------------------------

  interface str
    module procedure str_character
    module procedure str_integer
    module procedure str_real
    module procedure str_logical
    module procedure str_complex
  end interface
  
  interface char
    module procedure char_String
  end interface
  
  interface int
    module procedure int_String
  end interface
  
  interface dble
    module procedure dble_String
  end interface
  
  interface cmplx
    module procedure cmplx_String
  end interface
  
  interface lgcl
    module procedure lgcl_String
  end interface

  interface operator(//)
    module procedure concatenate_character_integer
    module procedure concatenate_integer_character
    module procedure concatenate_character_real
    module procedure concatenate_real_character
    module procedure concatenate_character_logical
    module procedure concatenate_logical_character
    module procedure concatenate_character_complex
    module procedure concatenate_complex_character
    
    module procedure concatenate_character_integers
    module procedure concatenate_integers_character
    module procedure concatenate_character_reals
    module procedure concatenate_reals_character
    module procedure concatenate_character_logicals
    module procedure concatenate_logicals_character
    module procedure concatenate_character_complexes
    module procedure concatenate_complexes_character
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
  
contains

! ----------------------------------------------------------------------
! Assignment
! ----------------------------------------------------------------------
! String = character(*)
pure subroutine assign_String_character(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  character(*),  intent(in)    :: that
  
  if(allocated(this%contents)) then
    deallocate(this%contents)
  endif
  
  allocate(character(len(that)) :: this%contents)
  this%contents = that
end subroutine

! String = String
pure subroutine assign_String_String(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  class(String), intent(in)    :: that
  
  ! An allocated() check is needed because casting
  !   from character(:) to character(*) failes if the character(:) is not
  !   allocated.
  if (allocated(that%contents)) then
    this = that%contents
  elseif (allocated(this%contents)) then
    deallocate(this%contents)
  endif
end subroutine

! String = integer
pure subroutine assign_String_integer(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  integer,       intent(in)    :: that
  
  character(12) :: temp
  
  write(temp,"(I0)") that
  
  this = trim(temp)
end subroutine

! String = real(dp)
pure subroutine assign_String_real(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  real(dp),      intent(in)    :: that
  
  integer, parameter :: width = 25
  integer, parameter :: decimal_places = 17
  type(String)       :: format_string
  
  character(width) :: temp
  
  format_string = str("(ES")//width//'.'//decimal_places//")"
  write(temp,char(format_string)) that
  
  this = temp
end subroutine

! String = logical
pure subroutine assign_String_logical(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  logical,       intent(in)    :: that
  
  if (that) then
    this = "T"
  else
    this = "F"
  endif
end subroutine

! String = logical
pure subroutine assign_String_complex(this,that)
  implicit none
  
  class(String), intent(inout) :: this
  complex(dp),   intent(in)    :: that
  
  type(String) :: imag
  
  imag = trim(str(aimag(that)))
  if (slice(imag,1,1)/='-') then
    imag = '+'//imag
  endif
  this = str(real(that))//imag//'i'
end subroutine

! character = String
pure subroutine assign_character_String(this,that)
  implicit none
  
  character(*),  intent(inout) :: this
  class(String), intent(in)    :: that
  
  this = that%contents
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

elemental function str_logical(this) result(output)
  implicit none
  
  logical, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

elemental function str_complex(this) result(output)
  implicit none
  
  complex(dp), intent(in) :: this
  type(String)            :: output
  
  output = this
end function

! ----------------------------------------------------------------------
! Conversion from String
! ----------------------------------------------------------------------
! character = char(String)
pure function char_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  character(len(this))     :: output
  
  output = this%contents
end function

! integer = int(String)
elemental function int_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  read(this%contents,*) output
end function

! real(dp) = dble(String)
elemental function dble_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  real(dp)                 :: output
  
  read(this%contents,*) output
end function

! complex(dp) = cmplx(String)
elemental function cmplx_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  complex(dp)              :: output
  
  integer :: i
  logical :: split_allowed
  
  split_allowed = .false.
  do i=len(this)-1,1,-1
    if (slice(this,i,i)=='E') then
      split_allowed = .true.
    elseif (split_allowed .and. slice(this,i,i)=='+') then
      output = cmplx(  dble(slice(this,1,i-1)), &
                    &  dble(slice(this,i+1,len(this)-1)), &
                    &  dp)
      exit
    elseif (split_allowed .and. slice(this,i,i)=='-') then
      output = cmplx(  dble(slice(this,1,i-1)), &
                    & -dble(slice(this,i+1,len(this)-1)), &
                    &  dp)
      exit
    endif
  enddo
end function

! logical = lgcl(String)
elemental function lgcl_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  logical                  :: output
  
  read(this%contents,*) output
end function

! ----------------------------------------------------------------------
! Concatenation
! ----------------------------------------------------------------------
! String = String//String
pure function concatenate_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = this%contents//that%contents
end function

! String = String//character
pure function concatenate_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  type(String)              :: output
  
  output = this%contents//that
end function

! String = character//String
pure function concatenate_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  type(String)              :: output
  
  output = this//that%contents
end function

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

! ----------------------------------------------------------------------
! Equality
! ----------------------------------------------------------------------
! String==String
elemental function equality_String_String(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = this%contents==that%contents
end function

! String==character
elemental function equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = this%contents==that
end function

! character==String
elemental function equality_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = this==that%contents
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
  
  output = this%contents/=that%contents
end function

! String/=character
elemental function non_equality_String_character(this,that) result(output)
  implicit none
  
  class(String), intent(in) :: this
  character(*),  intent(in) :: that
  logical                   :: output
  
  output = this%contents/=that
end function

! character/=String
elemental function non_equality_character_String(this,that) result(output)
  implicit none
  
  character(*),  intent(in) :: this
  class(String), intent(in) :: that
  logical                   :: output
  
  output = this/=that%contents
end function

! ----------------------------------------------------------------------
! Unary operators
! ----------------------------------------------------------------------
! integer = len(String)
elemental function len_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = len(this%contents)
end function

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
  integer      :: first     ! The position of a delimiter
  integer      :: second    ! The posiition of the next delimiter after 'first'
  integer      :: count     ! The number of tokens
  
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
  if (output%contents(1:1)/='-') then
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
pure function slice(this,first,last) result(output)
  implicit none
  
  type(String),  intent(in) :: this
  integer,       intent(in) :: first
  integer,       intent(in) :: last
  type(String)              :: output
  
  output = this%contents(first:last)
end function
end module
