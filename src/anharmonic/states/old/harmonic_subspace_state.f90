! ======================================================================
! A product of harmonic eigenstates along each mode in a degenerate subspace.
! ======================================================================
module harmonic_subspace_state_module
  use common_module
  
  use mode_state_module
  use subspace_state_module
  use harmonic_mode_basis_module
  use harmonic_mode_state_module
  implicit none
  
  private
  
  public :: HarmonicSubspaceState
  
  type, extends(Stringsable) :: HarmonicSubspaceState
    type(HarmonicModeState), allocatable :: states(:)
  contains
    ! I/O.
    procedure, public :: read  => read_HarmonicSubspaceState
    procedure, public :: write => write_HarmonicSubspaceState
  end type
  
  interface HarmonicSubspaceState
    module procedure new_HarmonicSubspaceState
    module procedure new_HarmonicSubspaceState_Strings
    module procedure new_HarmonicSubspaceState_StringArray
  end interface
  
  interface size
    module procedure size_HarmonicSubspaceState
  end interface
  
  !interface SubspaceState
  !  module procedure new_SubspaceState_HarmonicSubspaceState
  !end interface
contains

! Constructor and size() function.
function new_HarmonicSubspaceState(states) result(this)
  implicit none
  
  type(HarmonicModeState), intent(in) :: states(:)
  type(HarmonicSubspaceState)         :: this
  
  this%states = states
end function

function size_HarmonicSubspaceState(this) result(output)
  implicit none
  
  type(HarmonicSubspaceState), intent(in) :: this
  integer                                 :: output
  
  output = size(this%states)
end function

!! Convert a harmonic state to an explicitly represented state.
!function new_SubspaceState_HarmonicSubspaceState(this,bases) result(output)
!  implicit none
!  
!  type(HarmonicSubspaceState), intent(in) :: this
!  type(HarmonicModeBasis),     intent(in) :: bases(:)
!  type(SubspaceState)                     :: output
!  
!  type(ModeState), allocatable :: states(:)
!  type(HarmonicModeBasis)      :: basis
!  
!  integer :: i,ialloc
!  
!  allocate(states(size(this)), stat=ialloc); call err(ialloc)
!  do i=1,size(this)
!    basis = bases(first(bases%mode_id==this%states(i)%mode_id))
!    states(i) = ModeState(this%states(i), basis)
!  enddo
!  
!  output = SubspaceState(states)
!end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_HarmonicSubspaceState(this,input)
  implicit none
  
  class(HarmonicSubspaceState), intent(out) :: this
  type(String),                 intent(in)  :: input(:)
  
  select type(this); type is(HarmonicSubspaceState)
    this = HarmonicSubspaceState(HarmonicModeState(input))
  class default
    call err()
  end select
end subroutine

function write_HarmonicSubspaceState(this) result(output)
  implicit none
  
  class(HarmonicSubspaceState), intent(in) :: this
  type(String), allocatable                :: output(:)
  
  select type(this); type is(HarmonicSubspaceState)
    output = str(this%states)
  class default
    call err()
  end select
end function

function new_HarmonicSubspaceState_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)    :: input(:)
  type(HarmonicSubspaceState) :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicSubspaceState_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicSubspaceState)   :: this
  
  this = HarmonicSubspaceState(str(input))
end function
end module
