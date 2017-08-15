! ======================================================================
! Global constants.
! ======================================================================
module constants_module
  implicit none
  
  ! ----------------------------------------------------------------------
  ! Double precision.
  ! ----------------------------------------------------------------------
  integer, parameter :: dp = selected_real_kind(15,300)
  
  ! ----------------------------------------------------------------------
  ! Mathematical constants.
  ! ----------------------------------------------------------------------
  real(dp), parameter :: pi = 3.14159265358979324_dp
  
  ! ----------------------------------------------------------------------
  ! Cartesian directions.
  ! ----------------------------------------------------------------------
  character(1), parameter :: directions(3) = [ 'x','y','z' ]
  
  ! ----------------------------------------------------------------------
  ! Physical constants and unit conversions.
  ! ----------------------------------------------------------------------
  
  ! --------------------------------------------------
  ! Values taken from CODATA 2014, from NIST (physics.nist.gov) on 5/5/2017.
  ! --------------------------------------------------
  
  ! Hartree energy in eV.
  real(dp), parameter :: ev_per_hartree = 27.21138602_dp
  
  ! Boltzmann's constant in eV/K.
  real(dp), parameter :: kb_in_ev_per_k = 8.6173303e-5_dp
  
  ! Rydberg constant in eV.
  real(dp), parameter :: ev_per_rydberg = 13.605693009_dp
  
  ! Bohr radius in Angstrom.
  real(dp), parameter :: angstrom_per_bohr = 0.52917721067_dp
  
  ! Electron mass in kg.
  real(dp), parameter :: kg_per_me = 9.10938356e-31_dp
  
  ! Atomic mass unit in kg.
  real(dp), parameter :: kg_per_amu = 1.660539040e-27_dp
  
  ! Energy in cm-1 in eV.
  real(dp), parameter :: ev_per_inverse_cm = 2.1947463137e-4_dp
  
  ! --------------------------------------------------
  ! Derived values.
  ! --------------------------------------------------
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter :: kb_in_au = kb_in_ev_per_k / ev_per_rydberg
  
  ! ----------------------------------------------------------------------
  ! lte parameters.
  ! ----------------------------------------------------------------------
  ! Number of bins into which the frequency range is divided.
  integer, parameter :: max_bin = 1500
  
  ! Number of random samples of Brillouin zone to be made in each DoS set.
  integer, parameter :: no_samples = 50000
  
  ! Number of frequency DoS sets (for calculating error bars)
  integer, parameter :: no_fdos_sets = 20
  
end module
