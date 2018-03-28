! ======================================================================
! A coefficient times a product of Univariates.
! ======================================================================
module monomial_module
  use common_module
  
  use coupling_module
  use univariate_module
  use mode_displacement_module
  implicit none
  
  private
  
  public :: Monomial
  public :: size
  public :: operator(*)
  public :: operator(/)
  
  type, extends(Stringable) :: Monomial
    complex(dp)                   :: coefficient
    type(Univariate), allocatable :: modes(:)
  contains
    procedure :: evaluate   => evaluate_Monomial
    procedure :: derivative => derivative_Monomial
    procedure :: str        => str_Monomial
  end type
  
  interface Monomial
    module procedure new_Monomial
  end interface
  
  interface size
    module procedure size_Monomial
  end interface
  
  interface operator(*)
    module procedure multiply_Monomial_complex
    module procedure multiply_complex_Monomial
  end interface
  
  interface operator(/)
    module procedure divide_Monomial_complex
  end interface
contains

! Constructs a monomial from a mode coupling.
function new_Monomial(input) result(output)
  implicit none
  
  type(CoupledModes), intent(in) :: input
  type(Monomial)                 :: output
  
  integer, allocatable :: mode_ids(:)
  
  integer :: i,ialloc
  
  mode_ids = input%ids
  mode_ids = mode_ids(set(mode_ids))
  mode_ids = mode_ids(sort(mode_ids))
  output%coefficient = 1
  allocate(output%modes(size(mode_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(mode_ids)
    output%modes(i) = Univariate( id    = mode_ids(i), &
                                & power = count(input%ids==mode_ids(i)))
  enddo
end function

! The number of variables in the monomial.
function size_Monomial(this) result(output)
  implicit none
  
  type(Monomial), intent(in) :: this
  integer                    :: output
  
  output = size(this%modes)
end function

! Evaluates a Monomial at a given displacement.
function evaluate_Monomial(this,displacement) result(output)
  implicit none
  
  class(Monomial),        intent(in) :: this
  type(ModeDisplacement), intent(in) :: displacement
  complex(dp)                        :: output
  
  integer :: i,j
  
  output = this%coefficient
  
  do i=1,size(this)
    ! Find the mode in the displacement which matches that in the monomial.
    j = first(displacement%displacements%id==this%modes(i)%id)
    
    ! If the mode is not present in the displacement, then the displacement
    !    is zero. As such, the monomial is zero. (0**n=0 if n>0).
    if (j==0) then
      output = 0.0_dp
      return
    endif
    
    ! If the mode is present in both, evaluate the univariate at the
    !    displacement.
    output = output * this%modes(i)%evaluate(displacement%displacements(j))
  enddo
end function

! Takes the derivative of the Monomial in the direction of the given mode.
function derivative_Monomial(this,mode_id) result(output)
  implicit none
  
  class(Monomial), intent(in) :: this
  integer,         intent(in) :: mode_id
  type(Monomial)              :: output
  
  integer :: i
  
  ! Find the univariate corresponding to the mode.
  i = first(this%modes%id==mode_id)
  
  if (i==0) then
    ! If the mode is not present, then the derivative is zero.
    output%coefficient = 0
    output%modes = [Univariate::]
  else
    ! If the mode is present, then u^n -> n*u^(n-1).
    output = this
    output%coefficient = output%coefficient * output%modes(i)%power
    output%modes(i)%power = output%modes(i)%power - 1
    
    ! If n-1=0, remove that univariate.
    if (output%modes(i)%power==0) then
      output%modes = [output%modes(:i-1), output%modes(i+1:)]
    endif
  endif
end function

! Multiplication and division.
impure elemental function multiply_Monomial_complex(this,that) result(output)
  implicit none
  
  type(Monomial), intent(in) :: this
  complex(dp),    intent(in) :: that
  type(Monomial)             :: output
  
  output = this
  output%coefficient = output%coefficient * that
end function

impure elemental function multiply_complex_Monomial(this,that) result(output)
  implicit none
  
  complex(dp),    intent(in) :: this
  type(Monomial), intent(in) :: that
  type(Monomial)             :: output
  
  output = that
  output%coefficient = output%coefficient * this
end function

impure elemental function divide_Monomial_complex(this,that) result(output)
  implicit none
  
  type(Monomial), intent(in) :: this
  complex(dp),    intent(in) :: that
  type(Monomial)             :: output
  
  output = this
  output%coefficient = output%coefficient / that
end function

! I/O.
function str_Monomial(this) result(output)
  implicit none
  
  class(Monomial), intent(in) :: this
  type(String)                :: output
  
  integer :: i
  
  output = this%coefficient
  do i=1,size(this%modes)
    output = output//'.'//str(this%modes(i))
  enddo
end function
end module
