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
  
  ! Hartree energy in eV.
  real(dp), parameter, public :: ev_per_hartree = 27.21138602_dp
  
  ! Boltzmann's constant in eV/K.
  real(dp), parameter, public :: kb_in_ev_per_k = 8.6173303e-5_dp
  
  ! Rydberg constant in eV.
  real(dp), parameter, public :: ev_per_rydberg = 13.605693009_dp
  
  ! Bohr radius in Angstrom.
  real(dp), parameter, public :: angstrom_per_bohr = 0.52917721067_dp
  
  ! Electron mass in kg.
  real(dp), parameter, public :: kg_per_me = 9.10938356e-31_dp
  
  ! Atomic mass unit in kg.
  real(dp), parameter, public :: kg_per_amu = 1.660539040e-27_dp
  
  ! Energy in cm-1 in eV.
  real(dp), parameter, public :: ev_per_inverse_cm = 2.1947463137e-4_dp
  
  ! --------------------------------------------------
  ! Derived values.
  ! --------------------------------------------------
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter, public :: kb_in_au = kb_in_ev_per_k / ev_per_rydberg
end module
