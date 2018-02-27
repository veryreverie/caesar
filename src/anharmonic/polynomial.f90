! ======================================================================
! A sum of monomials.
! ======================================================================
module polynomial_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use monomial_module
  use stringable_module
  use printable_module
  implicit none
  
  private
  
  public :: Polynomial
  public :: size
  
  type, extends(Printable) :: Polynomial
    type(Monomial), allocatable :: terms(:)
  contains
    procedure :: evaluate   => evaluate_Polynomial
    procedure :: derivative => derivative_Polynomial
    procedure :: str        => str_Polynomial
  end type
  
  interface size
    module procedure size_Polynomial
  end interface
contains

! The number of terms in the polynomial.
function size_Polynomial(this) result(output)
  implicit none
  
  class(Polynomial), intent(in) :: this
  integer                       :: output
  
  output = size(this%terms)
end function

! Evaluates a Polynomial at a given displacement.
function evaluate_Polynomial(this,displacement) result(output)
  use mode_displacement_module
  implicit none
  
  class(Polynomial),      intent(in) :: this
  type(ModeDisplacement), intent(in) :: displacement
  complex(dp)                        :: output
  
  integer :: i
  
  output = 0
  
  do i=1,size(this)
    output = output + this%terms(i)%evaluate(displacement)
  enddo
end function

! Takes the derivative of the Polynomial in the direction of the given mode.
function derivative_Polynomial(this,mode_id) result(output)
  use logic_module
  implicit none
  
  class(Polynomial), intent(in) :: this
  integer,           intent(in) :: mode_id
  type(Polynomial)              :: output
  
  integer :: i,ialloc
  
  allocate(output%terms(size(this)), stat=ialloc); call err(ialloc)
  
  ! Take derivatives, term by term.
  do i=1,size(this)
    output%terms(i) = this%terms(i)%derivative(mode_id)
  enddo
  
  ! Remove the terms which are now zero.
  output%terms = output%terms(filter(output%terms,non_zero))
contains
  ! A Lambda for identifying non-zero monomials.
  function non_zero(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(Monomial)
      output = size(input)/=0
    end select
  end function
end function

! I/O.
function str_Polynomial(this) result(output)
  implicit none
  
  class(Polynomial), intent(in) :: this
  type(String), allocatable     :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  if (size(this)>=1) then
    output(1) = '  '//str(this%terms(1))
  endif
  do i=2,size(this)
    output(i) = '+ '//str(this%terms(i))
  enddo
end function
end module
