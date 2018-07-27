! ======================================================================
! A harmonic eigenstate along a single real normal mode.
! ======================================================================
module single_harmonic_state_module
  use common_module
  
  use single_mode_state_module
  use single_harmonic_basis_module
  implicit none
  
  private
  
  public :: SingleHarmonicState
  public :: SingleModeState
  
  ! The state |i> along mode j, with frequency w.
  type, extends(Stringable) :: SingleHarmonicState
    integer  :: mode_id    ! j.
    integer  :: occupation ! i.
  contains
    ! I/O.
    procedure, public :: read  => read_SingleHarmonicState
    procedure, public :: write => write_SingleHarmonicState
  end type
  
  interface SingleHarmonicState
    module procedure new_SingleHarmonicState
    module procedure new_SingleHarmonicState_String
  end interface
  
  interface SingleModeState
    module procedure new_SingleModeState_SingleHarmonicState
  end interface
contains

! Constructor.
function new_SingleHarmonicState(mode_id,occupation) result(this)
  implicit none
  
  integer, intent(in)       :: mode_id
  integer, intent(in)       :: occupation
  type(SingleHarmonicState) :: this
  
  this%mode_id    = mode_id
  this%occupation = occupation
end function

! Convert a harmonic state to a single mode state.
function new_SingleModeState_SingleHarmonicState(this,basis) result(output)
  implicit none
  
  type(SingleHarmonicState), intent(in) :: this
  type(SingleHarmonicBasis), intent(in) :: basis
  type(SingleModeState)                 :: output
  
  output = basis%state(this%occupation)
  
  if (output%mode_id/=this%mode_id) then
    call print_line(ERROR//': State modes do not match.')
    call err()
  endif
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_SingleHarmonicState(this,input)
  implicit none
  
  class(SingleHarmonicState), intent(out) :: this
  type(String),               intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  integer :: mode_id
  integer :: occupation
  
  select type(this); type is(SingleHarmonicState)
    line = split_line(input, delimiter='_')
    
    mode_id = int(slice(line(2),1,len(line(2))-1))
    occupation = int(slice(line(1),2,len(line(1))))
    
    this = SingleHarmonicState( mode_id    = mode_id,   &
                              & occupation = occupation )
  class default
    call err()
  end select
end subroutine

function write_SingleHarmonicState(this) result(output)
  implicit none
  
  class(SingleHarmonicState), intent(in) :: this
  type(String)                           :: output
  
  select type(this); type is(SingleHarmonicState)
    output = '|'//this%occupation//'_'//this%mode_id//'>'
  class default
    call err()
  end select
end function

impure elemental function new_SingleHarmonicState_String(input) result(this)
  implicit none
  
  type(String), intent(in)  :: input
  type(SingleHarmonicState) :: this
  
  this = input
end function
end module
