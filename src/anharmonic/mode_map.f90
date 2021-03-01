! ======================================================================
! A map of the potential along a given mode.
! ======================================================================
module caesar_mode_map_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: ModeMap
  
  type, extends(Stringsable) :: ModeMap
    integer               :: mode_id
    real(dp)              :: harmonic_frequency
    real(dp), allocatable :: mode_displacements(:)
    real(dp), allocatable :: l2_cartesian_displacements(:)
    real(dp), allocatable :: harmonic_energies(:)
    real(dp), allocatable :: harmonic_forces(:)
    real(dp), allocatable :: anharmonic_energies(:)
    real(dp), allocatable :: anharmonic_forces(:)
    real(dp), allocatable :: anharmonic_pressures(:)
    real(dp), allocatable :: sampled_energies(:)
    real(dp), allocatable :: sampled_forces(:)
    real(dp), allocatable :: sampled_pressures(:)
  contains
    procedure, public :: read  => read_ModeMap
    procedure, public :: write => write_ModeMap
  end type
  
  interface ModeMap
    module function new_ModeMap(mode_id,harmonic_frequency,               &
       & mode_displacements,l2_cartesian_displacements,harmonic_energies, &
       & harmonic_forces,anharmonic_energies,anharmonic_forces,           &
       & anharmonic_pressures,sampled_energies,sampled_forces,            &
       & sampled_pressures) result(this) 
      integer,  intent(in)           :: mode_id
      real(dp), intent(in)           :: harmonic_frequency
      real(dp), intent(in)           :: mode_displacements(:)
      real(dp), intent(in)           :: l2_cartesian_displacements(:)
      real(dp), intent(in)           :: harmonic_energies(:)
      real(dp), intent(in)           :: harmonic_forces(:)
      real(dp), intent(in), optional :: anharmonic_energies(:)
      real(dp), intent(in), optional :: anharmonic_forces(:)
      real(dp), intent(in), optional :: anharmonic_pressures(:)
      real(dp), intent(in), optional :: sampled_energies(:)
      real(dp), intent(in), optional :: sampled_forces(:)
      real(dp), intent(in), optional :: sampled_pressures(:)
      type(ModeMap)                  :: this
    end function
  
    module function new_ModeMap_potential(mode_displacements,              &
       & l2_cartesian_displacements,mode,potential,stress,anharmonic_data) &
       & result(this) 
      real(dp),             intent(in)           :: mode_displacements(:)
      real(dp),             intent(in)           :: l2_cartesian_displacements(:)
      type(RealMode),       intent(in)           :: mode
      class(PotentialData), intent(in), optional :: potential
      class(StressData),    intent(in), optional :: stress
      type(AnharmonicData), intent(in), optional :: anharmonic_data
      type(ModeMap)                              :: this
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module subroutine read_ModeMap(this,input) 
      class(ModeMap), intent(out) :: this
      type(String),   intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_ModeMap(this) result(output) 
      class(ModeMap), intent(in) :: this
      type(String), allocatable  :: output(:)
    end function
  end interface
  
  interface ModeMap
    module function new_ModeMap_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(ModeMap)            :: this
    end function
  
    impure elemental module function new_ModeMap_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(ModeMap)                 :: this
    end function
  end interface
end module
