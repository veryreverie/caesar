! ======================================================================
! A generalised harmonic basis function.
! ======================================================================
module basis_function_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  ! A one-mode term, e.g.
  !    ( a + b(u1) + c(u1)**2 )*e**(E*(u1)**2) => exponent=E, coeffs=[a,b,c]
  type, private :: OneModeTerm
    real(dp), allocatable :: exponent
    real(dp), allocatable :: coefficients(:)
  contains
    procedure :: evaluate => evaluate_OneModeTerm_real
  end type
  
  ! A product of one-mode terms.
  type :: BasisFunction
    type(OneModeTerm), allocatable, private :: terms(:)
  contains
    procedure :: evaluate => evaluate_BasisFunction_SamplingPoint
  end type
contains

! ------------------------------------------------------------
! Evaluates the one-mode term at a given displacement along the normal mode.
! ------------------------------------------------------------
function evaluate_OneModeTerm_real(this,u) result(output)
  implicit none
  
  class(OneModeTerm), intent(in) :: this
  real(dp),           intent(in) :: u
  real(dp)                       :: output
  
  integer :: i
  
  output = 0
  do i=1,size(this%coefficients)
    output = output + this%coefficients(i)*u**(i-1)
  enddo
  output = output * exp(this%exponent*u*u)
end function

! ------------------------------------------------------------
! Evaluates the basis function at a given sampling point.
! ------------------------------------------------------------
function evaluate_BasisFunction_SamplingPoint(this,sample_point) &
   & result(output)
  use sampling_point_module
  implicit none
  
  class(BasisFunction), intent(in) :: this
  type(SamplingPoint),  intent(in) :: sample_point
  real(dp)                         :: output
  
  integer :: i
  
  output = 1
  do i=1,size(this%terms)
    output = output * this%terms(i)%evaluate(sample_point%coordinates(i))
  enddo
end function

function generate_basis_functions() result(basis_functions)
  implicit none
  
  type(BasisFunction), allocatable :: basis_functions(:)
end function
end module
