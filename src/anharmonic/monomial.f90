! ======================================================================
! Monomials in terms of complex normal mode co-ordinates.
! ======================================================================
module monomial_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use stringable_module
  implicit none
  
  private
  
  ! A monomial, e.g.
  !    C * (u1)**a * (u2)**b * (u4)**d => coef=C, powers=[a,b,0,d]
  type, public, extends(Stringable) :: Monomial
    real(dp), public              :: coefficient
    integer,  public, allocatable :: powers(:)
  contains
    ! Evaluate energy or forces at a given value of u.
    procedure, public :: evaluate_energy
    procedure, public :: evaluate_forces
    
    ! I/O.
    procedure, public, pass(that) :: assign_String => assign_String_Monomial
  end type
contains

! ----------------------------------------------------------------------
! Evaluates the energy of a Monomial at a given displacement.
! ----------------------------------------------------------------------
function evaluate_energy(this,displacement) result(output)
  use normal_mode_module
  implicit none
  
  class(Monomial),  intent(in) :: this
  type(ModeVector), intent(in) :: displacement
  real(dp)                     :: output
  
  integer :: i
  
  output = this%coefficient
  do i=1,size(this%powers)
    output = output * displacement%vector(i)**this%powers(i)
  enddo
end function

! ----------------------------------------------------------------------
! Evaluates the force of a Monomial at a given displacement.
! Returns the result in normal mode co-ordinates.
! ----------------------------------------------------------------------
function evaluate_forces(this,displacement) result(output)
  use normal_mode_module
  implicit none
  
  class(Monomial),  intent(in) :: this
  type(ModeVector), intent(in) :: displacement
  type(ModeVector)             :: output
  
  integer        :: no_modes
  type(Monomial) :: derivative
  
  integer :: i,ialloc
  
  no_modes = size(displacement%vector)
  allocate(output%vector(no_modes), stat=ialloc); call err(ialloc)
  do i=1,no_modes
    derivative = this
    if (derivative%powers(i)==0) then
      derivative%coefficient = 0
      derivative%coefficient = 0
    else
      derivative%coefficient = derivative%coefficient &
                           & * derivative%powers(i)
      derivative%powers(i) = derivative%powers(i) - 1
    endif
    output%vector(i) = derivative%evaluate_energy(displacement)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine assign_String_Monomial(this,that)
  implicit none
  
  type(String),    intent(inout) :: this
  class(Monomial), intent(in)    :: that
  
  integer :: i
  
  this = that%coefficient
  do i=1,size(that%powers)
    if (that%powers(i)/=0) then
      this = this//'*u'//i//'^'//that%powers(i)
    endif
  enddo
end subroutine
end module
