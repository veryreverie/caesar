! ======================================================================
! A simple heap-allocated string class
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
  
  ! Binary operators
  public :: operator(//)  ! Concatenation between Strings and Strings/chars/ints
  
  ! Comparison operators
  public :: operator(==)  ! Equality comparison between Strings and Strings/chars
  
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
    module procedure assign_string_character
    module procedure assign_string_string
    module procedure assign_string_integer
    module procedure assign_character_string
  end interface

  interface drop
    module procedure drop_string
  end interface

  interface str
    module procedure str_character
    module procedure str_integer
  end interface

  interface char
    module procedure char_string
  end interface

  interface operator(//)
    module procedure concatenate_string_string
    module procedure concatenate_string_character
    module procedure concatenate_character_string
    module procedure concatenate_string_integer
    module procedure concatenate_integer_string
  end interface
  
  interface operator(==)
    module procedure equality_string_string
    module procedure equality_string_character
    module procedure equality_character_string
  end interface
  
  interface len
    module procedure len_string
  end interface
  
  interface system
    module procedure system_string
  end interface

contains

! ----------------------------------------------------------------------
! Assignment
! ----------------------------------------------------------------------
! String = character(*)
pure subroutine assign_string_character(output,input)
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
pure subroutine assign_string_string(output,input)
  implicit none
  
  type(String), intent(in)    :: input
  type(String), intent(inout) :: output
  
  output = input%contents
end subroutine

! String = integer
pure subroutine assign_string_integer(output,input)
  implicit none
  
  integer,      intent(in)    :: input
  type(String), intent(inout) :: output
  
  character(12) :: temp
  
  write(temp,"(I0)") input
  
  output = trim(temp)
end subroutine

! character = String
pure subroutine assign_character_string(output,input)
  implicit none
  
  type(String), intent(in)    :: input
  character(*), intent(inout) :: output
  
  output = input%contents
end subroutine

! ----------------------------------------------------------------------
! Deallocation
! ----------------------------------------------------------------------
pure subroutine drop_string(this)
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

! ----------------------------------------------------------------------
! Conversion from String
! ----------------------------------------------------------------------
! character = char(String)
pure function char_string(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  character(len(this))     :: output
  
  output = this%contents
end function

! ----------------------------------------------------------------------
! Concatenation
! ----------------------------------------------------------------------
! String = String//String
pure function concatenate_string_string(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b%contents
end function

! String = String//character
pure function concatenate_string_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  type(String)             :: output
  
  output = a%contents//b
end function

! String = character//String
pure function concatenate_character_string(a,b) result(output)
  implicit none
  
  character(*), intent(in) :: a
  type(String), intent(in) :: b
  type(String)             :: output
  
  output = a//b%contents
end function

! String = String//integer
pure function concatenate_string_integer(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  integer,      intent(in) :: b
  type(String) :: output
  
  type(String) :: temp
  
  temp = b
  
  output = a%contents//temp
end function

! String = String//integer
pure function concatenate_integer_string(a,b) result(output)
  implicit none
  
  integer,      intent(in) :: a
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
pure function equality_string_string(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  type(String), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b%contents
end function

! String==character
pure function equality_string_character(a,b) result(output)
  implicit none
  
  type(String), intent(in) :: a
  character(*), intent(in) :: b
  logical                  :: output
  
  output = a%contents==b
end function

! character==String
pure function equality_character_string(a,b) result(output)
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
pure function len_string(this) result(output)
  implicit none
  
  type(String), intent(in) :: this
  integer                  :: output
  
  output = len(this%contents)
end function

subroutine system_string(this)
  implicit none
  
  type(String), intent(in) :: this
  
  call system(char(this))
end subroutine

end module
