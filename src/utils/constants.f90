! ------------------------------------------------------------
! Global constants.
! ------------------------------------------------------------
module constants
  implicit none
  
  ! ----------------------------------------
  ! Double precision.
  ! ----------------------------------------
  integer, parameter :: dp = selected_real_kind(15,300)
  
  ! ----------------------------------------
  ! Identity matrix.
  ! ----------------------------------------
  integer, parameter :: identity(3,3) = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  
  ! ----------------------------------------
  ! Mathematical constants.
  ! ----------------------------------------
  real(dp), parameter :: pi = 3.14159265358979324_dp
  
  ! ----------------------------------------
  ! Cartesian directions.
  ! ----------------------------------------
  character(1), parameter :: directions(3) = (/ 'x','y','z' /)
  
  ! ----------------------------------------
  ! Physical constants and unit conversions.
  ! ----------------------------------------
  
  ! Values taken from CODATA 2014, from NIST (physics.nist.gov) on 5/5/2017.
  
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
  
  ! Derived values.
  
  ! Boltzmann's constant in Hartree per Kelvin.
  real(dp), parameter :: kb_in_au = kb_in_ev_per_k / ev_per_rydberg
  
  
! ! n.b. eV was defined variously as 27.211396132 and 27.21138602
! real(dp), parameter :: eV = 27.211396132_dp     ! Hartree energy eV
! real(dp), parameter :: thermal = 3.1577464e5_dp ! eV/kB
! real(dp), parameter :: kB = 8.6173324e-5_dp     ! Boltzmann's constant eV/K
! real(dp), parameter :: eV_per_A_to_au = 0.01944689814638725057_dp ! = bohr/eV
! real(dp), parameter :: kB_au_per_K = 3.16679002948702e-6_dp ! = kB/eV
! real(dp), parameter :: Ry = 13.605698066_dp     ! Rydberg energy eV
! real(dp), parameter :: bohr = 0.52911721092_dp  ! Bohr radius in Angstrom
!  
!  ! Parameters taken from Wikipedia. TODO: find a better source.
!  real(dp), parameter :: kg_to_me = 1.0977691e-32_dp ! kg to electron mass
!  real(dp), parameter :: amu_to_me = 1.822888486192e3_dp ! amu to electron mass
  
  ! ----------------------------------------
  ! lte parameters.
  ! ----------------------------------------
  ! Number of bins into which the frequency range is divided.
  integer, parameter :: max_bin = 1500
  
  ! Number of random samples of Brillouin zone to be made in each DoS set.
  integer, parameter :: no_samples = 50000
  
  ! Number of frequency DoS sets (for calculating error bars)
  integer, parameter :: no_fdos_sets = 20
  
end module
