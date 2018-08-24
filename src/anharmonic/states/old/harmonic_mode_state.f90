! ======================================================================
! A harmonic eigenstate along a single complex normal mode.
! ======================================================================
module harmonic_mode_state_module
  use common_module
  
  use mode_state_module
  use harmonic_mode_basis_module
  implicit none
  
  private
  
  public :: HarmonicModeState
  public :: ModeState
  
  ! The state |i> along mode j, with frequency w.
  type, extends(Stringable) :: HarmonicModeState
    integer  :: mode_id    ! j.
    integer  :: occupation ! i.
  contains
    procedure, public :: bloch_momentum
    ! I/O.
    procedure, public :: read  => read_HarmonicModeState
    procedure, public :: write => write_HarmonicModeState
  end type
  
  interface HarmonicModeState
    module procedure new_HarmonicModeState
    module procedure new_HarmonicModeState_String
  end interface
  
  interface ModeState
    module procedure new_ModeState_HarmonicModeState
  end interface
contains

! Constructor.
function new_HarmonicModeState(mode_id,occupation) result(this)
  implicit none
  
  integer, intent(in)     :: mode_id
  integer, intent(in)     :: occupation
  type(HarmonicModeState) :: this
  
  this%mode_id    = mode_id
  this%occupation = occupation
end function

! ----------------------------------------------------------------------
! Returns the Bloch momentum of the state.
! ----------------------------------------------------------------------
! The state |i> at q-point q has a Bloch momentum of i*q.
function bloch_momentum(this,mode,qpoints) result(output)
  implicit none
  
  class(HarmonicModeState), intent(in) :: this
  type(ComplexMode),        intent(in) :: mode
  type(QpointData),         intent(in) :: qpoints(:)
  type(FractionVector)                 :: output
  
  type(QpointData) :: qpoint
  
  if (this%mode_id/=mode%id) then
    call print_line(CODE_ERROR//': Trying to find the Bloch momentum of an &
       &incompatible state and mode.')
    call err()
  endif
  
  qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
  
  output = qpoint%qpoint * this%occupation
end function

! Convert a harmonic state to an explicitly represented state.
function new_ModeState_HarmonicModeState(this,basis) result(output)
  implicit none
  
  type(HarmonicModeState), intent(in) :: this
  type(HarmonicModeBasis), intent(in) :: basis
  type(ModeState)                     :: output
  
  if (basis%mode_id/=this%mode_id) then
    call print_line(ERROR//': Basis and state correspond to different modes.')
    call err()
  endif
  
  output = basis%state(this%occupation)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicModeState(this,input)
  implicit none
  
  class(HarmonicModeState), intent(out) :: this
  type(String),             intent(in)  :: input
  
  type(String), allocatable :: line(:)
  
  integer :: mode_id
  integer :: occupation
  
  select type(this); type is(HarmonicModeState)
    line = split_line(input, delimiter='_')
    
    mode_id = int(slice(line(2),1,len(line(2))-1))
    occupation = int(slice(line(1),2,len(line(1))))
    
    this = HarmonicModeState( mode_id    = mode_id,   &
                            & occupation = occupation )
  class default
    call err()
  end select
end subroutine

function write_HarmonicModeState(this) result(output)
  implicit none
  
  class(HarmonicModeState), intent(in) :: this
  type(String)                         :: output
  
  select type(this); type is(HarmonicModeState)
    output = '|'//this%occupation//'_'//this%mode_id//'>'
  class default
    call err()
  end select
end function

impure elemental function new_HarmonicModeState_String(input) result(this)
  implicit none
  
  type(String), intent(in) :: input
  type(HarmonicModeState)  :: this
  
  call this%read(input)
end function
end module
