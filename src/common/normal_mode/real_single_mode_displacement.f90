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
    procedure, public :: read  => read_RealSingleModeDisplacement
    procedure, public :: write => write_RealSingleModeDisplacement
  end type
  
  interface RealSingleModeDisplacement
    module procedure new_RealSingleModeDisplacement
  end interface
contains

! Constructor.
function new_RealSingleModeDisplacement(id,displacement) result(this)
  implicit none
  
  integer,  intent(in)             :: id
  real(dp), intent(in)             :: displacement
  type(RealSingleModeDisplacement) :: this
  
  this%id           = id
  this%displacement = displacement
end function

! I/O.
subroutine read_RealSingleModeDisplacement(this,input)
  implicit none
  
  class(RealSingleModeDisplacement), intent(out) :: this
  type(String),                      intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  real(dp)                  :: displacement
  
  select type(this); type is(RealSingleModeDisplacement)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse real single mode displacement &
         &from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1 then split_string = ["u3","=","2.1"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    displacement = dble(split_string(3))
    
    this = RealSingleModeDisplacement(id,displacement)
  end select
end subroutine

function write_RealSingleModeDisplacement(this) result(output)
  implicit none
  
  class(RealSingleModeDisplacement), intent(in) :: this
  type(String)                                  :: output
  
  select type(this); type is(RealSingleModeDisplacement)
    output = 'u'//this%id//' = '//this%displacement
  end select
end function
end module