! ======================================================================
! Harmonic states.
! ======================================================================
module harmonic_states_module
  use common_module
  
  use anharmonic_common_module
  implicit none
  
  private
  
  public :: startup_harmonic_states
  
  public :: HarmonicStates
  
  type, extends(BasisStates) :: HarmonicStates
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStates
    ! I/O.
    procedure, public :: read  => read_HarmonicStates
    procedure, public :: write => write_HarmonicStates
  end type
  
  interface HarmonicStates
    module procedure new_HarmonicStates
    module procedure new_HarmonicStates_Strings
    module procedure new_HarmonicStates_StringArray
  end interface
contains

! Startup procedure and type representation.
subroutine startup_harmonic_states
  implicit none
  
  type(HarmonicStates) :: states
  
  call states%startup()
end subroutine

impure elemental function representation_HarmonicStates() result(output)
  implicit none
  
  type(String) :: output
  
  output = 'harmonic state'
end function

! Constructor.
impure elemental function new_HarmonicStates() result(this)
  implicit none
  
  type(HarmonicStates) :: this
end function

! I/O.
subroutine read_HarmonicStates(this,input)
  implicit none
  
  class(HarmonicStates), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  select type(this); type is(HarmonicStates)
    this = HarmonicStates()
  class default
    call err()
  end select
end subroutine

function write_HarmonicStates(this) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  integer :: ialloc
  
  select type(this); type is(HarmonicStates)
    allocate(output(0), stat=ialloc); call err(ialloc)
  class default
    call err()
  end select
end function

function new_HarmonicStates_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(HarmonicStates)     :: this
  
  call this%read(input)
end function

impure elemental function new_HarmonicStates_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(HarmonicStates)          :: this
  
  this = HarmonicStates(str(input))
end function
end module
