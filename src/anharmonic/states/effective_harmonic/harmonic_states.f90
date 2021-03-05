! ======================================================================
! Harmonic states.
! ======================================================================
module caesar_harmonic_states_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
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
  
  interface
    impure elemental module function representation_HarmonicStates() &
       & result(output) 
      type(String) :: output
    end function
  end interface
  
  interface HarmonicStates
    ! Constructors.
    impure elemental module function new_HarmonicStates(subspace_id, &
       & frequency,thermal_energy) result(this) 
      integer,  intent(in) :: subspace_id
      real(dp), intent(in) :: frequency
      real(dp), intent(in) :: thermal_energy
      type(HarmonicStates) :: this
    end function
  
    recursive module function new_HarmonicStates_BasisStates(input) &
       & result(this) 
      class(BasisStates), intent(in) :: input
      type(HarmonicStates)           :: this
    end function
  end interface
  
  interface
    ! Cast a class(BasisStates) to a pointer of type(HarmonicStates).
    ! N.B. this must only be called on inputs with the TARGET attribute.
    recursive module function harmonic_states_pointer(input) result(this) 
      class(BasisStates), intent(in), target :: input
      type(HarmonicStates), pointer          :: this
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_HarmonicStates(this,input) 
      class(HarmonicStates), intent(out) :: this
      type(String),          intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_HarmonicStates(this) result(output) 
      class(HarmonicStates), intent(in) :: this
      type(String), allocatable         :: output(:)
    end function
  end interface
  
  interface HarmonicStates
    module function new_HarmonicStates_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(HarmonicStates)     :: this
    end function
  
    impure elemental module function new_HarmonicStates_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(HarmonicStates)          :: this
    end function
  end interface
end module
