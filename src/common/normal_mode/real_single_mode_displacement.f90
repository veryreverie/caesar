! ======================================================================
! A displacement along a single real mode.
! ======================================================================
module real_single_mode_displacement_submodule
  use utils_module
  implicit none
  
  private
  
  public :: RealSingleModeDisplacement
  
  ! The displacement along a single complex mode.
  type, extends(Stringable) :: RealSingleModeDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The displacement along the mode.
    real(dp) :: displacement
  contains
    procedure, public :: to_String => to_String_RealSingleModeDisplacement
  end type
  
  interface RealSingleModeDisplacement
    module procedure new_RealSingleModeDisplacement_String
  end interface
contains

! I/O.
function to_String_RealSingleModeDisplacement(this) result(output)
  implicit none
  
  class(RealSingleModeDisplacement), intent(in) :: this
  type(String)                                  :: output
  
  output = 'u'//this%id//' = '//this%displacement
end function

function new_RealSingleModeDisplacement_String(input) result(this)
  implicit none
  
  type(String), intent(in)         :: input
  type(RealSingleModeDisplacement) :: this
  
  type(String), allocatable :: split_string(:)
  
  split_string = split(input)
  if (size(split_string)/=3) then
    call print_line(ERROR//': unable to parse real single mode displacement &
       &from string: '//input)
    call err()
  endif
  
  ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
  ! The 'u' needs stripping off the first element to give the id.
  this%id = int(slice(split_string(1),2,len(split_string(1))))
  this%displacement = dble(split_string(3))
end function
end module
