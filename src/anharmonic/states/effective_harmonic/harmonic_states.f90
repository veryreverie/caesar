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
  
  public :: harmonic_states_pointer
  
  type, extends(BasisStates) :: HarmonicStates
    real(dp) :: frequency
    real(dp) :: thermal_energy
  contains
    procedure, public, nopass :: representation => &
                               & representation_HarmonicStates
    ! I/O.
    procedure, public :: read  => read_HarmonicStates
    procedure, public :: write => write_HarmonicStates
  end type
  
  interface HarmonicStates
    module procedure new_HarmonicStates
    module procedure new_HarmonicStates_BasisStates
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

! Constructors.
impure elemental function new_HarmonicStates(subspace_id,frequency, &
   & thermal_energy) result(this)
  implicit none
  
  integer,  intent(in) :: subspace_id
  real(dp), intent(in) :: frequency
  real(dp), intent(in) :: thermal_energy
  type(HarmonicStates) :: this
  
  this%subspace_id    = subspace_id
  this%frequency      = frequency
  this%thermal_energy = thermal_energy
end function

recursive function new_HarmonicStates_BasisStates(input) result(this)
  implicit none
  
  class(BasisStates), intent(in) :: input
  type(HarmonicStates)           :: this
  
  select type(input); type is(HarmonicStates)
    this = input
  type is(BasisStatesPointer)
    ! WORKAROUND: ifort doesn't recognise the interface to this function
    !    from within this function, so the full name is used instead.
    this = new_HarmonicStates_BasisStates(input%states())
  class default
    call err()
  end select
end function

! Cast a class(BasisStates) to a pointer of type(HarmonicStates).
! N.B. this must only be called on inputs with the TARGET attribute.
recursive function harmonic_states_pointer(input) result(this)
  implicit none
  
  class(BasisStates), intent(in), target :: input
  type(HarmonicStates), pointer          :: this
  
  select type(input); type is(HarmonicStates)
    this => input
  type is(BasisStatesPointer)
    this => harmonic_states_pointer(input%states_pointer())
  class default
    call err()
  end select
end function

! I/O.
subroutine read_HarmonicStates(this,input)
  implicit none
  
  class(HarmonicStates), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  integer  :: subspace_id
  real(dp) :: frequency
  real(dp) :: thermal_energy
  
  select type(this); type is(HarmonicStates)
    subspace_id = int(token(input(1),3))
    frequency = dble(token(input(2),3))
    thermal_energy = dble(token(input(3),4))
    
    this = HarmonicStates(subspace_id, frequency, thermal_energy)
  class default
    call err()
  end select
end subroutine

function write_HarmonicStates(this) result(output)
  implicit none
  
  class(HarmonicStates), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(HarmonicStates)
    output = [ 'Subspace       : '//this%subspace_id,   &
             & 'Frequency      : '//this%frequency,     &
             & 'Thermal energy : '//this%thermal_energy ]
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
