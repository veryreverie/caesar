! ======================================================================
! A state in a wavevector basis.
! See wavevector_basis.f90 for more information.
! ======================================================================
module wavevector_state_module
  use common_module
  implicit none
  
  private
  
  public :: WavevectorState
  
  type, extends(Stringable) :: WavevectorState
    real(dp), allocatable :: coefficients(:)
  contains
    ! I/O.
    procedure, public :: read  => read_WavevectorState
    procedure, public :: write => write_WavevectorState
  end type
  
  interface WavevectorState
    module procedure new_WavevectorState
    module procedure new_WavevectorState_String
  end interface
contains

! Constructor.
function new_WavevectorState(coefficients) result(this)
  implicit none
  
  real(dp), intent(in)  :: coefficients(:)
  type(WavevectorState) :: this
  
  this%coefficients = coefficients
end function

! I/O.
subroutine read_WavevectorState(this,input)
  implicit none
  
  class(WavevectorState), intent(out) :: this
  type(String),           intent(in)  :: input
  
  real(dp), allocatable :: coefficients(:)
  
  select type(this); type is(WavevectorState)
    coefficients = dble(split_line(input))
    
    this = WavevectorState(coefficients)
  end select
end subroutine

function write_WavevectorState(this) result(output)
  implicit none
  
  class(WavevectorState), intent(in) :: this
  type(String)                       :: output
  
  select type(this); type is(WavevectorState)
    output = join(this%coefficients, delimiter=' ')
  end select
end function

impure elemental function new_WavevectorState_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(WavevectorState)    :: this
  
  call this%read(input)
end function
end module
