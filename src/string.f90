! ======================================================================
! A simple heap-allocated String class
! ======================================================================
module string_module
  implicit none
  
  private
  
  ! ----------------------------------------------------------------------
  ! Public interface
  ! ----------------------------------------------------------------------
  
  ! String class
  public :: String
  
  ! Allocate and deallocate
  public :: assignment(=) ! Assignment to and from String
  public :: drop          ! Deallocator
  
  ! Conversions between classes
  public :: str           ! Conversion to String
  public :: char          ! Conversion from String to character
  public :: int           ! Conversion from String to integer
  public :: dble          ! Conversion from String to real(dp)
  
  ! Binary operators
  public :: operator(//)
  
  ! Comparison operators
  public :: operator(==)
  public :: operator(/=)
  
  ! Unary operators
  public :: len           ! character-like len
  public :: lower_case    ! convert to lower case
  public :: split         ! split by spaces
  
  ! Operators with side-effects
  public :: system
  
  type String
    character(:), allocatable, private :: contents
  end type
  
  ! ----------------------------------------------------------------------
  ! Interfaces
  ! ----------------------------------------------------------------------

  interface assignment(=)
    module procedure assign_String_character
    module procedure assign_String_String
    module procedure assign_String_integer
    module procedure assign_String_real
    module procedure assign_String_logical
    module procedure assign_character_String
  end interface

  interface drop
    module procedure drop_String
  end interface

  interface str
    module procedure str_character
    module procedure str_integer
    module procedure str_real
    module procedure str_logical
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

  interface operator(//)
    module procedure concatenate_String_String
    module procedure concatenate_String_character
    module procedure concatenate_character_String
    module procedure concatenate_String_integer
    module procedure concatenate_integer_String
    module procedure concatenate_String_real
    module procedure concatenate_real_String
    module procedure concatenate_String_logical
    module procedure concatenate_logical_String
  end interface
  
  interface operator(==)
    module procedure equality_String_String
    module procedure equality_String_character
    module procedure equality_character_String
  end interface
  
  interface operator(/=)
    module procedure non_equality_String_String
    module procedure non_equality_String_character
    module procedure non_equality_character_String
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
  
  interface system
    module procedure system_String
  end interface

contains

! ----------------------------------------------------------------------
! Assignment
! ----------------------------------------------------------------------
! String = character(*)
pure subroutine assign_String_character(output,input)
  implicit none
  
  character(*), intent(in)    :: input
  type(String), intent(inout) :: output
  
  if(allocated(output%contents)) then
    deallocate(output%contents)
  endif
  
  allocate(character(len(input)) :: output%contents)
  output%contents = input
end subroutine

! String = String
pure subroutine assign_String_String(output,input)
  implicit none
  
  type(String), intent(in)    :: input
  type(String), intent(inout) :: output
  
  output = input%contents
end subroutine

! String = integer
pure subroutine assign_String_integer(output,input)
  implicit none
  
  integer,      intent(in)    :: input
  type(String), intent(inout) :: output
  
  character(12) :: temp
  
  write(temp,"(I0)") input
  
  output = trim(temp)
end subroutine

! String = real(dp)
pure subroutine assign_String_real(output,input)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in)    :: input
  type(String), intent(inout) :: output
  
  integer, parameter :: width = 25
  integer, parameter :: decimal_places = 17
  type(String)       :: format_string
  
  character(width) :: temp
  
  format_string = str("(ES")//width//'.'//decimal_places//")"
  write(temp,char(format_string)) input
  
  output = temp
end subroutine

! String = logical
pure subroutine assign_String_logical(output,input)
  implicit none
  
  logical,      intent(in)    :: input
  type(String), intent(inout) :: output
  
  if (input) then
    output = "T"
  else
    output = "F"
  endif
end subroutine

! character = String
pure subroutine assign_character_String(output,input)
  implicit none
  
  type(String), intent(in)    :: input
  character(*), intent(inout) :: output
  
  output = input%contents
end subroutine

! ----------------------------------------------------------------------
! Deallocation
! ----------------------------------------------------------------------
pure subroutine drop_String(this)
  implicit none
  
  type(String), intent(inout) :: this
  
  deallocate(this%contents)
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
  use constants, only : dp
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
  use constants, only : dp
  implicit none
  
  type(String), intent(in) :: this
  real(dp)                 :: output
  
  read(this%contents,*) output
end function

! ----------------------------------------------------------------------
! Concatenation
! ----------------------------------------------------------------------
! String = String//String
elemental function concatenate_String_String(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b%contents
end function

! String = String//character
elemental function concatenate_String_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b
end function

! String = character//String
elemental function concatenate_character_String(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a//b%contents
end function

! String = String//integer
elemental function concatenate_String_integer(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  integer,      intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = b
  
  output = a%contents//temp
end function

! String = integer//String
elemental function concatenate_integer_String(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = a
  
  output = temp//b%contents
end function

! String = String//real(dp)
elemental function concatenate_String_real(a,b) result(output)
  use constants, only : dp
  implicit none
  
  type(String), intent(in) :: a
  real(dp),     intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = b
  
  output = a%contents//temp
end function

! String = real(dp)//String
elemental function concatenate_real_String(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = a
  
  output = temp//b%contents
end function

! String = String//logical
elemental function concatenate_String_logical(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  logical,      intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = b
  
  output = a%contents//temp
end function

! String = logical//String
elemental function concatenate_logical_String(a,b) result(output)
  implicit none
  
  logical,      intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = a
  
  output = temp//b%contents
end function

! ----------------------------------------------------------------------
! Equality
! ----------------------------------------------------------------------
! String==String
elemental function equality_String_String(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b%contents
end function

! String==character
elemental function equality_String_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b
end function

! character==String
elemental function equality_character_String(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a==b%contents
end function

! ----------------------------------------------------------------------
! Non-equality
! ----------------------------------------------------------------------
! String==String
elemental function non_equality_String_String(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a%contents/=b%contents
end function

! String==character
elemental function non_equality_String_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  logical                  :: output
  
  output = a%contents/=b
end function

! character==String
elemental function non_equality_character_String(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a/=b%contents
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

! call system(String)
subroutine system_String(this)
  implicit none
  
  type(String), intent(in) :: this
  
  call system(char(this))
end subroutine

end module
