! ------------------------------------------------------------
! global constants
! ------------------------------------------------------------
module constants
  implicit none
  
  ! double precision
  integer,  parameter :: dp = kind(1.d0) ! defines double precision
  
  ! mathematical constants
  real(dp), parameter :: pi = 3.14159265358979324d0
  real(dp), parameter :: twopi = 2.d0*pi
  real(dp), parameter :: third = 1.d0/3.d0
  
  ! physical constants and unit conversions
  ! n.b. eV was defined variously as 27.211396132 and 27.21138602
  real(dp), parameter :: eV = 27.211396132d0   ! Hartree energy eV
  real(dp), parameter :: thermal = 3.1577464E5 ! Hartree temperature K
  real(dp), parameter :: kB = 8.6173324E-5     ! Boltzmann's constant eV/K
  real(dp), parameter :: eV_per_A_to_au = 0.01944689814638725057d0

  
contains

end module
