! ======================================================================
! ======================================================================
module full_subspace_states_module
  use common_module
  
  use anharmonic_common_module
  
  use wavevector_state_module
  implicit none
  
  private
  
  public :: startup_full_subspace_states
  
  public :: FullSubspaceStates
  
  type, extends(BasisStates) :: FullSubspaceStates
    type(WavevectorState), allocatable :: states(:)
    real(dp),              allocatable :: energies(:)
  contains
    procedure, public, nopass :: representation => &
                               & representation_FullSubspaceStates
    
    ! I/O.
    procedure, public :: read  => read_FullSubspaceStates
    procedure, public :: write => write_FullSubspaceStates
  end type
  
  interface FullSubspaceStates
    module procedure new_FullSubspaceStates
    module procedure new_FullSubspaceStates_BasisStates
    module procedure new_FullSubspaceStates_Strings
    module procedure new_FullSubspaceStates_StringArray
  end interface
contains

! Startup procedure.
subroutine startup_full_subspace_states()
  implicit none
  
  type(FullSubspaceStates) :: states
  
  call states%startup()
end subroutine

! ----------------------------------------------------------------------
! FullSubspaceStates methods.
! ----------------------------------------------------------------------
! Constructor.
function new_FullSubspaceStates(states,energies) result(this)
  implicit none
  
  type(WavevectorState), intent(in) :: states(:)
  real(dp),              intent(in) :: energies(:)
  type(FullSubspaceStates)          :: this
  
  this%states = states
  this%energies = energies
end function

recursive function new_FullSubspaceStates_BasisStates(input) result(this)
  implicit none
  
  class(BasisStates), intent(in) :: input
  type(FullSubspaceStates)       :: this
  
  select type(input); type is(FullSubspaceStates)
    this = input
  type is(BasisStatesPointer)
    this = FullSubspaceStates(input%states())
  class default
    call err()
  end select
end function

! Type representation.
impure elemental function representation_FullSubspaceStates() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'full_subspace'
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_FullSubspaceStates(this,input)
  implicit none
  
  class(FullSubspaceStates), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  type(StringArray), allocatable :: sections(:)
  
  type(WavevectorState), allocatable :: states(:)
  real(dp),              allocatable :: energies(:)
  
  type(String), allocatable :: line(:)
  
  integer :: i,ialloc
  
  select type(this); type is(FullSubspaceStates)
    sections = split_into_sections(input, separating_line=repeat('=',50))
    
    allocate( states(size(sections)),   &
            & energies(size(sections)), &
            & stat=ialloc); call err(ialloc)
    do i=1,size(sections)
      states(i) = WavevectorState(sections(i)%strings(:size(sections(i))-1))
      
      line = split_line(sections(i)%strings(size(sections(i))))
      energies(i) = dble(line(3))
    enddo
    this = FullSubspaceStates(states, energies)
  class default
    call err()
  end select
end subroutine

function write_FullSubspaceStates(this) result(output)
  implicit none
  
  class(FullSubspaceStates), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  type(StringArray), allocatable :: sections(:)
  
  integer :: i
  
  select type(this); type is(FullSubspaceStates)
    sections = [( StringArray([ str(this%states(i)),             &
                &               'Energy : '//this%energies(i)]), &
                & i=1,                                           &
                & size(this%states)                              )]
    output = str(sections, separating_line=repeat('=',50))
  class default
    call err()
  end select
end function

function new_FullSubspaceStates_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(FullSubspaceStates) :: this
  
  call this%read(input)
end function

impure elemental function new_FullSubspaceStates_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(FullSubspaceStates)      :: this
  
  this = FullSubspaceStates(str(input))
end function
end module
