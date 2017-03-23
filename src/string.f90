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
  public :: len        ! Character-like len.
  public :: lower_case ! Convert to lower case.
  public :: split      ! Split into String(:) by spaces
  public :: join       ! Join into single String by spaces.
  public :: pad_str    ! left pads integers without '-' signs with a ' '
  
  ! Operators with side-effects
  public :: system              ! Makes a system call.
  public :: read_line_from_user ! Reads a line from the terminal.
  public :: print_line          ! write(*,'(a)')
  
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
    
    module procedure concatenate_character_integer
    module procedure concatenate_integer_character
    module procedure concatenate_character_real
    module procedure concatenate_real_character
    module procedure concatenate_character_logical
    module procedure concatenate_logical_character
    
    module procedure concatenate_String_integers
    module procedure concatenate_integers_String
    module procedure concatenate_String_reals
    module procedure concatenate_reals_String
    module procedure concatenate_String_logicals
    module procedure concatenate_logicals_String
    
    module procedure concatenate_character_integers
    module procedure concatenate_integers_character
    module procedure concatenate_character_reals
    module procedure concatenate_reals_character
    module procedure concatenate_character_logicals
    module procedure concatenate_logicals_character
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
  
  interface join
    module procedure join_String
    module procedure join_real
    module procedure join_integer
    module procedure join_logical
  end interface
  
  interface system
    module procedure system_String
  end interface
  
  interface print_line
    module procedure print_line_character
    module procedure print_line_file_character
    module procedure print_line_String
    module procedure print_line_file_String
    
    module procedure print_line_integer
    module procedure print_line_file_integer
    module procedure print_line_real
    module procedure print_line_file_real
    module procedure print_line_logical
    module procedure print_line_file_logical
    
    module procedure print_line_integers
    module procedure print_line_file_integers
    module procedure print_line_reals
    module procedure print_line_file_reals
    module procedure print_line_logicals
    module procedure print_line_file_logicals
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
pure function concatenate_String_String(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b%contents
end function

! String = String//character
pure function concatenate_String_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b
end function

! String = character//String
pure function concatenate_character_String(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a//b%contents
end function

! String = String//integer
pure function concatenate_String_integer(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  integer,      intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = integer//String
pure function concatenate_integer_String(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = String//real(dp)
pure function concatenate_String_real(a,b) result(output)
  use constants, only : dp
  implicit none
  
  type(String), intent(in) :: a
  real(dp),     intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = real(dp)//String
pure function concatenate_real_String(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = String//logical
pure function concatenate_String_logical(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  logical,      intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = logical//String
pure function concatenate_logical_String(a,b) result(output)
  implicit none
  
  logical,      intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = character//integer
pure function concatenate_character_integer(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  integer,      intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = integer//character
pure function concatenate_integer_character(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a
  character(*), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = character//real(dp)
pure function concatenate_character_real(a,b) result(output)
  use constants, only : dp
  implicit none
  
  character(*), intent(in) :: a
  real(dp),     intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = real(dp)//character
pure function concatenate_real_character(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a
  character(*), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = character//logical
pure function concatenate_character_logical(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  logical,      intent(in) :: b
  type(String) :: output
  
  output = a//str(b)
end function

! String = logical//character
pure function concatenate_logical_character(a,b) result(output)
  implicit none
  
  logical,      intent(in) :: a
  character(*), intent(in) :: b
  type(String) :: output
  
  output = str(a)//b
end function

! String = String//integer(:)
pure function concatenate_String_integers(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  integer,      intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = integer(:)//String
pure function concatenate_integers_String(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a(:)
  type(String), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
end function

! String = String//real(dp)(:)
pure function concatenate_String_reals(a,b) result(output)
  use constants, only : dp
  implicit none
  
  type(String), intent(in) :: a
  real(dp),     intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = real(dp)(:)//String
pure function concatenate_reals_String(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a(:)
  type(String), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
end function

! String = String//logical(:)
pure function concatenate_String_logicals(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  logical,      intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = logical(:)//String
pure function concatenate_logicals_String(a,b) result(output)
  implicit none
  
  logical,      intent(in) :: a(:)
  type(String), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
end function

! String = character//integer(:)
pure function concatenate_character_integers(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  integer,      intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = integer(:)//character
pure function concatenate_integers_character(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a(:)
  character(*), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
end function

! String = character//real(dp)(:)
pure function concatenate_character_reals(a,b) result(output)
  use constants, only : dp
  implicit none
  
  character(*), intent(in) :: a
  real(dp),     intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = real(dp)(:)//character
pure function concatenate_reals_character(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a(:)
  character(*), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
end function

! String = character//logical(:)
pure function concatenate_character_logicals(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  logical,      intent(in) :: b(:)
  type(String) :: output
  
  output = a//join(b)
end function

! String = logical(:)//character
pure function concatenate_logicals_character(a,b) result(output)
  implicit none
  
  logical,      intent(in) :: a(:)
  character(*), intent(in) :: b
  type(String) :: output
  
  output = join(a)//b
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

! Joins a String(:) array into one String.
pure function join_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this(:)
  type(String)             :: output
  
  integer :: i
  
  if (size(this)==0) then
    output = ''
  else
    output = this(1)
    do i=2,size(this)
      output = output//' '//this(i)
    enddo
  endif
end function

! Converts real(:) to String(:) and then joins.
pure function join_real(this) result(output)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: this(:)
  type(String)         :: output
  
  output = join(str(this))
end function

! Converts integer(:) to String(:), pads positive numbers with a space, and joins.
pure function join_integer(this) result(output)
  implicit none
  
  integer, intent(in) :: this(:)
  type(String)        :: output
  
  output = join(pad_str(this))
end function

! Converts logical(:) to String(:) and then joins.
pure function join_logical(this) result(output)
  implicit none
  
  logical, intent(in) :: this(:)
  type(String)        :: output
  
  output = join(str(this))
end function

elemental function pad_str(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
  if (output%contents(1:1)/='-') then
    output = ' '//output
  endif
end function

! call system(String)
subroutine system_String(this)
  implicit none
  
  type(String), intent(in) :: this
  
  call system(char(this))
end subroutine

function read_line_from_user() result(line)
  implicit none
  
  type(String) :: line
  
  character(1000) :: char_line
  
  read(*,'(a)') char_line
  line = trim(char_line)
end function

! ----------------------------------------------------------------------
! Subroutines to print lines, as write(), but with error checking
!    and formatting.
! ----------------------------------------------------------------------
subroutine print_line_character(line)
  use err_module
  implicit none
  
  character(*), intent(in) :: line
  
  integer :: ierr
  
  write(*,'(a)',iostat=ierr) line
  if (ierr /= 0) then
    write(*,*) 'Error in print_line.'
    call err()
  endif
end subroutine

subroutine print_line_file_character(file_unit,line)
  use err_module
  implicit none
  
  integer,      intent(in) :: file_unit
  character(*), intent(in) :: line
  
  integer :: ierr
  
  write(file_unit,'(a)',iostat=ierr) line
  if (ierr /= 0) then
    write(*,*) 'Error in print_line.'
    call err()
  endif
end subroutine

subroutine print_line_String(line)
  implicit none
  
  type(String), intent(in) :: line
  
  call print_line(char(line))
end subroutine

subroutine print_line_file_String(file_unit,line)
  implicit none
  
  integer,      intent(in) :: file_unit
  type(String), intent(in) :: line
  
  call print_line(file_unit,char(line))
end subroutine

subroutine print_line_integer(this)
  implicit none
  
  integer, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_integer(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  integer, intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_real(this)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_real(file_unit,this)
  use constants, only : dp
  implicit none
  
  integer,  intent(in) :: file_unit
  real(dp), intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_logical(this)
  implicit none
  
  logical, intent(in) :: this
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_logical(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  logical, intent(in) :: this
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_integers(this)
  implicit none
  
  integer, intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_integers(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  integer, intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_reals(this)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_reals(file_unit,this)
  use constants, only : dp
  implicit none
  
  integer,  intent(in) :: file_unit
  real(dp), intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
end subroutine

subroutine print_line_logicals(this)
  implicit none
  
  logical, intent(in) :: this(:)
  
  call print_line(''//this)
end subroutine

subroutine print_line_file_logicals(file_unit,this)
  implicit none
  
  integer, intent(in) :: file_unit
  logical, intent(in) :: this(:)
  
  call print_line(file_unit,''//this)
end subroutine
end module
