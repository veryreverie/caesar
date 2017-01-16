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
  public :: operator(//)  ! Concatenation
  
  ! Comparison operators
  public :: operator(==)  ! Equality comparison
  
  ! Unary operators
  public :: len           ! character-like len
  
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
    module procedure assign_character_String
  end interface

  interface drop
    module procedure drop_String
  end interface

  interface str
    module procedure str_character
    module procedure str_integer
    module procedure str_real
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
  end interface
  
  interface operator(==)
    module procedure equality_String_String
    module procedure equality_String_character
    module procedure equality_character_String
  end interface
  
  interface len
    module procedure len_String
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
pure function str_character(this) result(output)
  implicit none
  
  character(*), intent(in) :: this
  type(String)             :: output
  
  output = this
end function

pure function str_integer(this) result(output)
  implicit none
  
  integer, intent(in) :: this
  type(String)        :: output
  
  output = this
end function

pure function str_real(this) result(output)
  use constants, only : dp
  implicit none
  
  real(dp), intent(in) :: this
  type(String)         :: output
  
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
pure function int_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  read(this%contents,*) output
end function

! real(dp) = dble(String)
pure function dble_String(this) result(output)
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
  
  type(String) :: temp
  
  temp = b
  
  output = a%contents//temp
end function

! String = integer//String
pure function concatenate_integer_String(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a
  type(String), intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = a
  
  output = temp//b%contents
end function

! String = String//real(dp)
pure function concatenate_String_real(a,b) result(output)
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
pure function concatenate_real_String(a,b) result(output)
  use constants, only : dp
  implicit none
  
  real(dp),     intent(in) :: a
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
pure function equality_String_String(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b%contents
end function

! String==character
pure function equality_String_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b
end function

! character==String
pure function equality_character_String(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a==b%contents
end function

! ----------------------------------------------------------------------
! Unary operators
! ----------------------------------------------------------------------
! integer = len(String)
pure function len_String(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = len(this%contents)
end function

subroutine system_String(this)
  implicit none
  
  type(String), intent(in) :: this
  
  call system(char(this))
end subroutine

end module