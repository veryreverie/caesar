! ======================================================================
! A base class for the String class.
! ======================================================================
! This class exists to provide an interface between the actual contents of a 
!    String and all of the procedures involving a String.
! This ensures that all of the memory management is in one place, and that the
!    String class cannot cause segfaults or allocation errors.
module string_base_submodule
  implicit none
  
  private
  
  public :: StringBase
  public :: char
  
  type :: StringBase
    character(:), allocatable, private :: contents_
  contains
    generic,   public  :: assignment(=) => assign_StringBase_character
    procedure, private ::                  assign_StringBase_character
  end type
  
  interface char
    module procedure char_StringBase
  end interface
contains

! ----------------------------------------------------------------------
! Public procedures.
! ----------------------------------------------------------------------

! Create a string from a character(*). Effectively a setter.
! Currently uses naive memory allocation, reallocating all memory for every
!    operation.
! If it is desired that String be sped up, this should probably be changed to
!    a smarter memory scheme, e.g. C++'s vector model.
pure subroutine assign_StringBase_character(output,input)
  implicit none
  
  class(StringBase), intent(inout) :: output
  character(*),      intent(in)    :: input
  
  output%contents_ = input
end subroutine

! Conversion to character(*). Effectively a getter.
function char_StringBase(this) result(output)
  use error_submodule
  implicit none
  
  class(StringBase), intent(in) :: this
  character(:), allocatable     :: output
  
  if (allocated(this%contents_)) then
    output = this%contents_
  else
    write(*,*) CODE_ERROR//': Trying to use the contents of a string before &
       &it has been allocated.'
    call abort_with_stacktrace()
  endif
end function
end module
