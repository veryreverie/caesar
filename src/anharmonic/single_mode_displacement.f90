! ======================================================================
! A displacement along a single complex mode.
! ======================================================================
module single_mode_displacement_module
  use common_module
  implicit none
  
  private
  
  public :: SingleModeDisplacement
  
  ! The displacement along a single complex mode.
  type, extends(Stringable) :: SingleModeDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The displacement along the mode.
    complex(dp) :: displacement
  contains
    procedure :: str => str_SingleModeDisplacement
  end type
contains

! I/O.
function str_SingleModeDisplacement(this) result(output)
  implicit none
  
  class(SingleModeDisplacement), intent(in) :: this
  type(String)                              :: output
  
  output = 'u'//this%id//' = '//this%displacement
end function
end module
