! ======================================================================
! A vector in cartesian co-ordinates.
! ======================================================================
module cartesian_vector_submodule
  use utils_module
  
  use structure_submodule
  implicit none
  
  private
  
  public :: CartesianVector
  public :: size
  public :: operator(*)
  public :: operator(+)
  public :: sum
  
  type, extends(Stringsable) :: CartesianVector
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_CartesianVector
    procedure, public :: write => write_CartesianVector
  end type
  
  interface CartesianVector
    module procedure new_CartesianVector
    module procedure new_CartesianVector_StringArray
  end interface
  
  interface size
    module procedure size_CartesianVector
  end interface
  
  interface operator(*)
    module procedure multiply_real_CartesianVector
    module procedure multiply_CartesianVector_real
  end interface
  
  interface operator(+)
    module procedure add_CartesianVector_CartesianVector
  end interface
  
  interface sum
    module procedure sum_CartesianVectors
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality: constructor and size() function.
! ----------------------------------------------------------------------
function new_CartesianVector(vectors) result(this)
  implicit none
  
  type(RealVector), intent(in) :: vectors(:)
  type(CartesianVector)        :: this
  
  this%vectors = vectors
end function

function size_CartesianVector(this) result(output)
  implicit none
  
  class(CartesianVector), intent(in) :: this
  integer                            :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
function multiply_real_CartesianVector(this,that) result(output)
  implicit none
  
  real(dp),              intent(in) :: this
  type(CartesianVector), intent(in) :: that
  type(CartesianVector)             :: output
  
  output = CartesianVector(this * that%vectors)
end function

function multiply_CartesianVector_real(this,that) result(output)
  implicit none
  
  type(CartesianVector), intent(in) :: this
  real(dp),              intent(in) :: that
  type(CartesianVector)             :: output
  
  output = CartesianVector(this%vectors * that)
end function

function add_CartesianVector_CartesianVector(this,that) &
   & result(output)
  implicit none
  
  type(CartesianVector), intent(in) :: this
  type(CartesianVector), intent(in) :: that
  type(CartesianVector)             :: output
  
  output = CartesianVector(this%vectors + that%vectors)
end function

function sum_CartesianVectors(input) result(output)
  implicit none
  
  type(CartesianVector), intent(in) :: input(:)
  type(CartesianVector)             :: output
  
  integer :: i
  
  if (size(input)==0) then
    call print_line(ERROR//': Trying to sum() an empty array.')
    call err()
  endif
  
  output = input(1)
  do i=2,size(input)
    output = output + input(i)
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianVector(this,input)
  implicit none
  
  class(CartesianVector), intent(out) :: this
  type(String),           intent(in)  :: input(:)
  
  select type(this); type is(CartesianVector)
    this = CartesianVector(RealVector(input))
  end select
end subroutine

function write_CartesianVector(this) result(output)
  implicit none
  
  class(CartesianVector), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  select type(this); type is(CartesianVector)
    output = str(this%vectors)
  end select
end function

impure elemental function new_CartesianVector_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianVector)         :: this
  
  this = input
end function
end module
