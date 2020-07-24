! ======================================================================
! A base class for the String class.
! ======================================================================
! This class exists to provide an interface between the actual contents of a 
!    String and all of the procedures involving a String.
! This ensures that all of the memory management is in one place, and that the
!    String class cannot cause segfaults or allocation errors.
module string_base_module
  use error_module
  implicit none
  
  private
  
  public :: StringBase
  public :: assignment(=)
  public :: char
  public :: operator(==)
  public :: operator(/=)
  
  type :: StringBase
    character(:), allocatable, private :: contents_
  contains
    procedure, private :: check
  end type
  
  interface assignment(=)
    module procedure assign_StringBase_character
  end interface
  
  interface char
    module procedure char_StringBase
  end interface
  
  interface operator(==)
    module procedure equality_StringBase_character
    module procedure equality_character_StringBase
    module procedure equality_StringBase_StringBase
  end interface
  
  interface operator(/=)
    module procedure non_equality_StringBase_character
    module procedure non_equality_character_StringBase
    module procedure non_equality_StringBase_StringBase
  end interface
contains

! --------------------------------------------------
! Checking that this StringBase is allocated before use.
! --------------------------------------------------
subroutine check(this)
  implicit none
  
  class(StringBase), intent(in) :: this
  
  if (.not. allocated(this%contents_)) then
    write(*,*) CODE_ERROR//': Trying to use the contents of a string before &
       &it has been allocated.'
    call err()
  endif
end subroutine

! --------------------------------------------------
! Assignment from character(*) and conversion to character(*).
! --------------------------------------------------

! Create a string from a character(*). Effectively a setter.
! Currently uses naive memory allocation, reallocating all memory for every
!    operation.
! If it is desired that String be sped up, this should probably be changed to
!    a smarter memory scheme, e.g. C++'s vector model.
subroutine assign_StringBase_character(output,input)
  implicit none
  
  class(StringBase), intent(out) :: output
  character(*),      intent(in)  :: input
  
  output%contents_ = input
end subroutine

! Conversion to character(*). Effectively a getter.
function char_StringBase(this) result(output)
  implicit none
  
  class(StringBase), intent(in) :: this
  character(:), allocatable     :: output
  
  call this%check()
  output = this%contents_
end function

! --------------------------------------------------
! Comparison with character(*) and StringBase.
! --------------------------------------------------
function equality_StringBase_character(this,that) result(output)
  implicit none
  
  class(StringBase), intent(in) :: this
  character(*),      intent(in) :: that
  logical                       :: output
  
  call this%check()
  output = this%contents_==that
end function

function equality_character_StringBase(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  class(StringBase), intent(in) :: that
  logical                       :: output
  
  call that%check()
  output = this==that%contents_
end function

impure elemental function equality_StringBase_StringBase(this,that) &
   & result(output)
  implicit none
  
  class(StringBase), intent(in) :: this
  class(StringBase), intent(in) :: that
  logical                       :: output
  
  call this%check()
  call that%check()
  output = this%contents_==that%contents_
end function

function non_equality_StringBase_character(this,that) result(output)
  implicit none
  
  class(StringBase), intent(in) :: this
  character(*),      intent(in) :: that
  logical                       :: output
  
  output = .not. this==that
end function

function non_equality_character_StringBase(this,that) result(output)
  implicit none
  
  character(*),      intent(in) :: this
  class(StringBase), intent(in) :: that
  logical                       :: output
  
  output = .not. this==that
end function

impure elemental function non_equality_StringBase_StringBase(this,that) &
   & result(output)
  implicit none
  
  class(StringBase), intent(in) :: this
  class(StringBase), intent(in) :: that
  logical                       :: output
  
  output = .not. this==that
end function
end module
