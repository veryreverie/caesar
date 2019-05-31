! ======================================================================
! ======================================================================
module wavevector_states_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_state_module
  implicit none
  
  private
  
  public :: startup_wavevector_states
  
  public :: WavevectorStates
  
  type, extends(BasisStates) :: WavevectorStates
    type(WavevectorState), allocatable :: states(:)
    real(dp),              allocatable :: energies(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_WavevectorStates
    
    ! I/O.
    procedure, public :: read  => read_WavevectorStates
    procedure, public :: write => write_WavevectorStates
  end type
  
  interface WavevectorStates
    module procedure new_WavevectorStates
    module procedure new_WavevectorStates_BasisStates
    module procedure new_WavevectorStates_Strings
    module procedure new_WavevectorStates_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_wavevector_states()
  implicit none
  
  type(WavevectorStates) :: states
  
  call states%startup()
end subroutine

! ----------------------------------------------------------------------
! WavevectorStates methods.
! ----------------------------------------------------------------------
! Constructor.
function new_WavevectorStates(states,energies) result(this)
  implicit none
  
  type(WavevectorState), intent(in) :: states(:)
  real(dp),              intent(in) :: energies(:)
  type(WavevectorStates)            :: this
  
  this%states = states
  this%energies = energies
end function

recursive function new_WavevectorStates_BasisStates(input) result(this)
  implicit none
  
  class(BasisStates), intent(in) :: input
  type(WavevectorStates)         :: this
  
  select type(input); type is(WavevectorStates)
    this = input
  type is(BasisStatesPointer)
    this = WavevectorStates(input%states())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_WavevectorStates() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'wavevector state'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_WavevectorStates(this,input)
  implicit none
  
  class(WavevectorStates), intent(out) :: this
  type(String),            intent(in)  :: input(:)
  
  type(StringArray), allocatable :: sections(:)
  
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(WavevectorStates)
    sections = split_into_sections(input, separating_line=repeat('=',50))
    
    allocate( states(size(sections)),   &
            & energies(size(sections)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(sections)
      states(i) = WavevectorState(sections(i)%strings(:size(sections(i))-1))
      
      line = split_line(sections(i)%strings(size(sections(i))))
      energies(i) = dble(line(3))
    enddo
    this = WavevectorStates(states, energies)
  class default
    call err()
  end select
end subroutine

function write_WavevectorStates(this) result(output)
  implicit none
  
  class(WavevectorStates), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i
  
  select type(this); type is(WavevectorStates)
    sections = [( StringArray([ str(this%states(i)),             &
                &               'Energy : '//this%energies(i)]), &
                & i=1,                                           &
                & size(this%states)                              )]
    output = str(sections, separating_line=repeat('=',50))
  class default
    call err()
  end select
end function

function new_WavevectorStates_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(WavevectorStates)   :: this
  
  call this%read(input)
end function

impure elemental function new_WavevectorStates_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(WavevectorStates)        :: this
  
  this = WavevectorStates(str(input))
end function
end module
