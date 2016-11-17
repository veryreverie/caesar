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
  
  ! physical constants
  real(dp), parameter :: eV = 27.211396132d0   ! Hartree energy / e
  real(dp), parameter :: thermal = 3.1577464E5 ! Hartree energy / kB
  
contains

end module
