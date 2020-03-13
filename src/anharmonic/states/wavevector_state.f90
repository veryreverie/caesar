! ======================================================================
! A state in a wavevector basis.
! See wavevector_basis.f90 for more information.
! ======================================================================
module wavevector_state_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: startup_wavevector_state
  
  public :: WavevectorState
  
  public :: wavevector_state_pointer
  
  type, extends(BasisState) :: WavevectorState
    type(FractionVector)  :: wavevector
    integer,  allocatable :: state_ids(:)
    real(dp), allocatable :: coefficients(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorState
    ! I/O.
    procedure, public :: read  => read_WavevectorState
    procedure, public :: write => write_WavevectorState
  end type
  
  interface WavevectorState
    module procedure new_WavevectorState
    module procedure new_WavevectorState_BasisState
    module procedure new_WavevectorState_Strings
    module procedure new_WavevectorState_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_wavevector_state()
  implicit none
  
  type(WavevectorState) :: state
  
  call state%startup()
end subroutine

! Type representation.
impure elemental function representation_WavevectorState() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'wavevector state'
end function

! Constructors.
function new_WavevectorState(subspace_id,wavevector,state_ids,coefficients) &
   & result(this) 
  implicit none
  
  integer,              intent(in) :: subspace_id
  type(FractionVector), intent(in) :: wavevector
  integer,              intent(in) :: state_ids(:)
  real(dp),             intent(in) :: coefficients(:)
  type(WavevectorState)            :: this
  
  this%subspace_id  = subspace_id
  this%wavevector   = wavevector
  this%state_ids    = state_ids
  this%coefficients = coefficients
end function

recursive function new_WavevectorState_BasisState(input) result(this)
  implicit none
  
  class(BasisState), intent(in) :: input
  type(WavevectorState)         :: this
  
  select type(input); type is(WavevectorState)
    this = input
  type is(BasisStatePointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    this = new_WavevectorState_BasisState(input%state())
  class default
    call err()
  end select
end function

! Cast a class(BasisState) to a pointer of type(WavevectorState).
! N.B. this must only be called on inputs with the TARGET attribute.
recursive function wavevector_state_pointer(input) result(this)
  implicit none
  
  class(BasisState), intent(in), target :: input
  type(WavevectorState), pointer        :: this
  
  select type(input); type is(WavevectorState)
    this => input
  type is(BasisStatePointer)
    this => wavevector_state_pointer(input%state())
  class default
    call err()
  end select
end function

! I/O.
subroutine read_WavevectorState(this,input)
  implicit none
  
  class(WavevectorState), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  integer               :: subspace_id
  type(FractionVector)  :: wavevector
  integer,  allocatable :: state_ids(:)
  real(dp), allocatable :: coefficients(:)
  
  select type(this); type is(WavevectorState)
    subspace_id = int(token(input(1), 3))
    wavevector = FractionVector(join(tokens(input(2),3)))
    state_ids = int(tokens(input(3),3))
    coefficients = dble(tokens(input(4),3))
    
    this = WavevectorState(subspace_id,wavevector,state_ids,coefficients)
  end select
end subroutine

function write_WavevectorState(this) result(output)
  implicit none
  
  class(WavevectorState), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(WavevectorState)
    output = [ 'Subspace     : '//this%subspace_id,                      &
             & 'Wavevector   : '//str(this%wavevector),                  &
             & 'States       : '//this%state_ids,                        &
             & 'Coefficients : '//join(this%coefficients, delimiter=' ') ]
  end select
end function

function new_WavevectorState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(WavevectorState)    :: this
  
  call this%read(input)
end function

impure elemental function new_WavevectorState_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(WavevectorState)         :: this
  
  this = WavevectorState(str(input))
end function
end module
