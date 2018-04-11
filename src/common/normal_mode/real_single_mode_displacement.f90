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
contains

! I/O.
function to_String_RealSingleModeDisplacement(this) result(output)
  implicit none
  
  class(RealSingleModeDisplacement), intent(in) :: this
  type(String)                                  :: output
  
  output = 'u'//this%id//' = '//this%displacement
end function
end module
