! ======================================================================
! A complex phase, of the form exp(2*pi*i*frac), where frac is an integer
!    fraction.
! ======================================================================
module phase_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use fraction_module
  implicit none
  
  type :: PhaseData
    type(IntFraction) :: fraction
  end type
  
  interface PhaseData
    module procedure new_PhaseData
  end interface
  
  interface cmplx
    module procedure cmplx_PhaseData
  end interface
contains

! Constructor.
function new_PhaseData(input) result(this)
  implicit none
  
  type(IntFraction), intent(in) :: input
  type(PhaseData)               :: this
  
  this%fraction = input
end function

! Conversion to complex(dp).
function cmplx_PhaseData(this) result(output)
  use constants_module, only : pi
  implicit none
  
  type(PhaseData) :: this
  complex(dp)     :: output
  
  real(dp) :: exponent
  
  exponent = 2*pi*dble(this%fraction)
  output = cmplx(cos(exponent),sin(exponent),dp)
end function

! Finds the exact phase of a complex number.
function calculate_phase(input,denom) result(output)
  use constants_module, only : pi
  implicit none
  
  complex(dp), intent(in) :: input
  integer,     intent(in) :: denom
  type(PhaseData)         :: output
  
  real(dp)          :: phase_real
  type(IntFraction) :: phase_frac
  
  phase_real = atan2(aimag(input), real(input)) / (2*pi)
  phase_frac = IntFraction(nint(phase_real*denom),denom)
  if (abs(dble(phase_frac)-phase_real)>0.01_dp) then
    call print_line(ERROR//': Phase incompatible with given denominator.')
    call err()
  endif
  output = PhaseData(phase_frac)
end function
end module
