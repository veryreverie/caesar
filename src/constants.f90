! ------------------------------------------------------------
! global constants
! ------------------------------------------------------------
module constants
  implicit none
  
  integer, parameter  :: dp = kind(1.d0) ! defines double precision
  real(dp), parameter :: pi = 3.14159265358979324d0
  real(dp), parameter :: twopi = 2.d0*pi
  real(dp), parameter :: third = 1.d0/3.d0
end module
