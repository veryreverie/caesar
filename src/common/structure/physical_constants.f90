! ======================================================================
! Global physical constants.
! ======================================================================
module physical_constants_module
  use utils_module
  implicit none
  
  private
  
  ! --------------------------------------------------
  ! Values taken from CODATA 2014, from NIST (physics.nist.gov) on 5/5/2017.
  ! --------------------------------------------------
  
  ! Twice Rydberg's constant in eV.
  real(dp), parameter, public :: EV_PER_HARTREE = 27.21138602_dp
  
  ! Boltzmann's constant in eV/K.
  real(dp), parameter, public :: KB_IN_EV_PER_K = 8.6173303e-5_dp
  
  ! Rydberg constant in eV.
  real(dp), parameter, public :: EV_PER_RYDBERG = 13.605693009_dp
  
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
  ! Derived values.
  ! --------------------------------------------------
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter, public :: KB_IN_AU = KB_IN_EV_PER_K / EV_PER_HARTREE
end module
