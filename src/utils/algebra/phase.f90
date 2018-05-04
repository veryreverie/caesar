! ======================================================================
! A complex phase, of the form exp(2*pi*i*frac), where frac is an integer
!    fraction.
! ======================================================================
module phase_submodule
  use precision_module
  use io_module
  
  use fraction_submodule
  use mathematical_constants_submodule
  implicit none
  
  private
  
  public :: PhaseData
  public :: cmplx
  public :: operator(==)
  public :: operator(/=)
  public :: calculate_phase
  
  type, extends(Stringable) :: PhaseData
    type(IntFraction) :: fraction
  contains
    procedure, public :: read  => read_PhaseData
    procedure, public :: write => write_PhaseData
  end type
  
  interface PhaseData
    module procedure new_PhaseData
  end interface
  
  interface cmplx
    module procedure cmplx_PhaseData
  end interface
  
  interface operator(==)
    module procedure equality_PhaseData_PhaseData
  end interface
  
  interface operator(/=)
    module procedure non_equality_PhaseData_PhaseData
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
  implicit none
  
  type(PhaseData) :: this
  complex(dp)     :: output
  
  real(dp) :: exponent
  
  exponent = 2*PI*dble(this%fraction)
  output = cmplx(cos(exponent),sin(exponent),dp)
end function

! Finds the exact phase of a complex number.
function calculate_phase(input,denom) result(output)
  implicit none
  
  complex(dp), intent(in) :: input
  integer,     intent(in) :: denom
  type(PhaseData)         :: output
  
  real(dp)          :: phase_real
  type(IntFraction) :: phase_frac
  
  phase_real = atan2(aimag(input), real(input)) / (2*PI)
  phase_frac = IntFraction(nint(phase_real*denom),denom)
  if (abs(dble(phase_frac)-phase_real)>0.01_dp) then
    call print_line(ERROR//': Phase incompatible with given denominator.')
    call print_line('Phase: '//phase_real)
    call print_line('Denominator: '//denom)
    call err()
  endif
  output = PhaseData(phase_frac)
end function

! ----------------------------------------------------------------------
! Comparison.
! ----------------------------------------------------------------------
impure elemental function equality_PhaseData_PhaseData(this,that) &
   & result(output)
  implicit none
  
  type(PhaseData), intent(in) :: this
  type(PhaseData), intent(in) :: that
  logical                     :: output
  
  output = this%fraction==that%fraction
end function

impure elemental function non_equality_PhaseData_PhaseData(this,that) &
   & result(output)
  implicit none
  
  type(PhaseData), intent(in) :: this
  type(PhaseData), intent(in) :: that
  logical                     :: output
  
  output = .not. this==that
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_PhaseData(this,input)
  implicit none
  
  class(PhaseData), intent(out) :: this
  type(String),     intent(in)  :: input
  
  select type(this); type is(PhaseData)
    this = PhaseData(frac(slice(input,10,len(input)-1)))
  end select
end subroutine

function write_PhaseData(this) result(output)
  implicit none
  
  class(PhaseData), intent(in) :: this
  type(String)                 :: output
  
  select type(this); type is(PhaseData)
    output = 'exp(2pii*'//this%fraction//')'
  end select
end function
end module
