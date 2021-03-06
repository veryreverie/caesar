! ======================================================================
! Global physical constants.
! ======================================================================
module caesar_physical_constants_module
  use caesar_utils_module
  implicit none
  
  private
  
  ! --------------------------------------------------
  ! Values taken from CODATA 2014, from NIST (physics.nist.gov) on 5/5/2017.
  ! --------------------------------------------------
  
  ! Twice Rydberg's constant in eV.
  real(dp), parameter, public :: EV_PER_HARTREE = 27.21138602_dp
  
  ! Boltzmann's constant in eV/K.
  real(dp), parameter, public :: KB_IN_EV_PER_K = 8.6173303e-5_dp
  
  ! Twice Rydberg's constant in J.
  real(dp), parameter, public :: JOULES_PER_HARTREE = 4.359744650e-18_dp
  
  ! Twice Rydberg's constant in inverse cm.
  real(dp), parameter, public :: INVERSE_CM_PER_HARTREE = 2.194746313702e4_dp
  
  ! Bohr radius in Angstrom.
  real(dp), parameter, public :: ANGSTROM_PER_BOHR = 0.52917721067_dp
  
  ! Electron mass in kg.
  real(dp), parameter, public :: KG_PER_ME = 9.10938356E-31_dp
  
  ! Atomic mass unit in kg.
  real(dp), parameter, public :: KG_PER_AMU = 1.660539040E-27_dp
  
  ! --------------------------------------------------
  ! Fundamental constants.
  ! --------------------------------------------------
  
  ! Rydberg's constant in Hartree.
  real(dp), parameter, public :: RYDBERG_PER_HARTREE = 2.0_dp
  
  real(dp), parameter, public :: RYDBERG_MASS_PER_ME = 0.5_dp
  
  ! --------------------------------------------------
  ! Derived values.
  ! --------------------------------------------------
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter, public :: KB_IN_AU = KB_IN_EV_PER_K / EV_PER_HARTREE
  
  ! AMU in ME.
  real(dp), parameter, public :: AMU_PER_ME = KG_PER_ME / KG_PER_AMU
  
  ! --------------------------------------------------
  ! Atomic symbols.
  ! --------------------------------------------------
  public :: atomic_symbol_to_number
  public :: atomic_number_to_symbol
  
  interface
    module subroutine initialise_atomic_symbols() 
    end subroutine
  end interface
  
  interface
    impure elemental module function atomic_number_to_symbol(input) &
       & result(output) 
      integer, intent(in) :: input
      type(String)        :: output
    end function
  end interface
  
  interface
    impure elemental module function atomic_symbol_to_number(input) &
       & result(output) 
      type(String), intent(in) :: input
      integer                  :: output
    end function
  end interface
end module
