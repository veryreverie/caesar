! ======================================================================
! A displacement along a single complex mode.
! ======================================================================
module complex_single_mode_displacement_submodule
  use utils_module
  implicit none
  
  private
  
  public :: ComplexSingleModeDisplacement
  
  ! The displacement along a single complex mode.
  type, extends(Stringable) :: ComplexSingleModeDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The displacement along the mode.
    complex(dp) :: displacement
  contains
    procedure, public :: to_String => to_String_ComplexSingleModeDisplacement
  end type
  
  interface ComplexSingleModeDisplacement
    module procedure new_ComplexSingleModeDisplacement_String
  end interface
contains

! I/O.
function to_String_ComplexSingleModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexSingleModeDisplacement), intent(in) :: this
  type(String)                                     :: output
  
  output = 'u'//this%id//' = '//this%displacement
end function

function new_ComplexSingleModeDisplacement_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(ComplexSingleModeDisplacement) ::this
  
  type(String), allocatable :: split_string(:)
  
  split_string = split(input)
  if (size(split_string)/=3) then
    call print_line(ERROR//': unable to parse complex single mode &
       &displacement from string: '//input)
    call err()
  endif
  
  ! If e.g. id=3 and power=2.1+1.2i then split_string = ["u3","=","2.1+1.2i"]
  ! The 'u' needs stripping off the first element to give the id.
  this%id = int(slice(split_string(1),2,len(split_string(1))))
  this%displacement = cmplx(split_string(3))
end function
end module
