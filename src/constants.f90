! ------------------------------------------------------------
! global constants
! ------------------------------------------------------------
module constants
  implicit none
  
  ! ----------------------------------------
  ! double precision
  ! ----------------------------------------
  integer, parameter :: dp = selected_real_kind(15,300)
  
  ! ----------------------------------------
  ! Identity matrix
  ! ----------------------------------------
  integer, parameter :: identity(3,3) = reshape((/1,0,0,0,1,0,0,0,1/),(/3,3/))
  
  ! ----------------------------------------
  ! mathematical constants
  ! ----------------------------------------
  real(dp), parameter :: pi = 3.14159265358979324_dp
  
  ! ----------------------------------------
  ! physical constants and unit conversions
  ! ----------------------------------------
  ! n.b. eV was defined variously as 27.211396132 and 27.21138602
  real(dp), parameter :: eV = 27.211396132_dp     ! Hartree energy eV
  real(dp), parameter :: thermal = 3.1577464e5_dp ! Hartree temperature K
  real(dp), parameter :: kB = 8.6173324e-5_dp     ! Boltzmann's constant eV/K
  real(dp), parameter :: eV_per_A_to_au = 0.01944689814638725057_dp
  real(dp), parameter :: kB_au_per_K = 3.16679002948702e-6_dp ! kB in Hartrees/K
  real(dp), parameter :: Ry = 13.605698066_dp     ! Rydberg energy eV
  real(dp), parameter :: bohr = 0.52911721092_dp  ! Bohr radius in Angstrom
  
  ! ----------------------------------------
  ! lte parameters
  ! ----------------------------------------
  ! Number of bins into which the frequency range is divided.
  integer, parameter :: max_bin = 1500
  
  ! Number of random samples of Brillouin zone to be made in each DoS set.
  integer, parameter :: no_samples = 50000
  
  ! Number of frequency DoS sets (for calculating error bars)
  integer, parameter :: no_fdos_sets = 20
  
end module
