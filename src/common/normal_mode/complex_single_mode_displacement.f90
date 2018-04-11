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
contains

! I/O.
function to_String_ComplexSingleModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexSingleModeDisplacement), intent(in) :: this
  type(String)                                     :: output
  
  output = 'u'//this%id//' = '//this%displacement
end function
end module
