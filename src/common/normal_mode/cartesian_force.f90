! ======================================================================
! A force in cartesian co-ordinates.
! ======================================================================
module cartesian_force_module
  use utils_module
  implicit none
  
  private
  
  public :: CartesianForce
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  
  type, extends(Stringsable) :: CartesianForce
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_CartesianForce
    procedure, public :: write => write_CartesianForce
  end type
  
  interface CartesianForce
    module procedure new_CartesianForce
    module procedure new_CartesianForce_Strings
    module procedure new_CartesianForce_StringArray
  end interface
  
  interface size
    module procedure size_CartesianForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_CartesianForce
    module procedure multiply_CartesianForce_real
  end interface
  
  interface operator(/)
    module procedure divide_CartesianForce_real
  end interface
  
  interface operator(+)
    module procedure add_CartesianForce_CartesianForce
  end interface
  
  interface sum
    module procedure sum_CartesianForces
  end interface
  
  interface operator(-)
    module procedure negative_CartesianForce
    module procedure subtract_CartesianForce_CartesianForce
  end interface
contains

! ----------------------------------------------------------------------
! Constructor and size() function.
! ----------------------------------------------------------------------
function new_CartesianForce(forces) result(this)
  implicit none
  
  type(RealVector), intent(in) :: forces(:)
  type(CartesianForce)         :: this
  
  this%vectors = forces
end function

function size_CartesianForce(this) result(output)
  implicit none
  
  class(CartesianForce), intent(in) :: this
  integer                           :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
impure elemental function multiply_real_CartesianForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),             intent(in) :: this
  type(CartesianForce), intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this * that%vectors)
end function

impure elemental function multiply_CartesianForce_real(this,that) &
   & result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  real(dp),             intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%vectors * that)
end function

impure elemental function divide_CartesianForce_real(this,that) result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  real(dp),             intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%vectors / that)
end function

impure elemental function add_CartesianForce_CartesianForce(this,that) &
   & result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  type(CartesianForce), intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%vectors + that%vectors)
end function

function sum_CartesianForces(input) result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: input(:)
  type(CartesianForce)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(ERROR//': Trying to sum an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

impure elemental function negative_CartesianForce(this) result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  type(CartesianForce)             :: output
  
  output = CartesianForce(-this%vectors)
end function

impure elemental function subtract_CartesianForce_CartesianForce(this,that) &
   & result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  type(CartesianForce), intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%vectors - that%vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianForce(this,input)
  implicit none
  
  class(CartesianForce), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  select type(this); type is(CartesianForce)
    this = CartesianForce(RealVector(input))
  class default
    call err()
  end select
end subroutine

function write_CartesianForce(this) result(output)
  implicit none
  
  class(CartesianForce), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(CartesianForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_CartesianForce_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(CartesianForce)     :: this
  
  call this%read(input)
end function

impure elemental function new_CartesianForce_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianForce)          :: this
  
  this = CartesianForce(str(input))
end function
end module
