! ======================================================================
! A vector in mass-weighted cartesian co-ordinates.
! ======================================================================
module mass_weighted_vector_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: MassWeightedVector
  public :: size
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: MassWeightedVector
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_MassWeightedVector
    procedure, public :: write => write_MassWeightedVector
  end type
  
  interface MassWeightedVector
    module procedure new_MassWeightedVector
    module procedure new_MassWeightedVector_StringArray
  end interface
  
  interface size
    module procedure size_MassWeightedVector
  end interface
  
  interface operator(*)
    module procedure multiply_real_MassWeightedVector
    module procedure multiply_MassWeightedVector_real
  end interface
  
  interface operator(/)
    module procedure divide_MassWeightedVector_real
  end interface
  
  interface operator(+)
    module procedure add_MassWeightedVector_MassWeightedVector
  end interface
  
  interface sum
    module procedure sum_MassWeightedVectors
  end interface
  
  interface operator(-)
    module procedure negative_MassWeightedVector
    module procedure subtract_MassWeightedVector_MassWeightedVector
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality: constructor and size() function.
! ----------------------------------------------------------------------
function new_MassWeightedVector(vectors) result(this)
  implicit none
  
  type(RealVector), intent(in) :: vectors(:)
  type(MassWeightedVector)     :: this
  
  this%vectors = vectors
end function

function size_MassWeightedVector(this) result(output)
  implicit none
  
  class(MassWeightedVector), intent(in) :: this
  integer                               :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
impure elemental function multiply_real_MassWeightedVector(this,that) &
   & result(output)
  implicit none
  
  real(dp),                 intent(in) :: this
  type(MassWeightedVector), intent(in) :: that
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(this * that%vectors)
end function

impure elemental function multiply_MassWeightedVector_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: this
  real(dp),                 intent(in) :: that
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(this%vectors * that)
end function

impure elemental function divide_MassWeightedVector_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: this
  real(dp),                 intent(in) :: that
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(this%vectors / that)
end function

impure elemental function add_MassWeightedVector_MassWeightedVector(this, &
   & that) result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: this
  type(MassWeightedVector), intent(in) :: that
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(this%vectors + that%vectors)
end function

function sum_MassWeightedVectors(input) result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: input(:)
  type(MassWeightedVector)             :: output
  
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

impure elemental function negative_MassWeightedVector(this) result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: this
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(-this%vectors)
end function

impure elemental function subtract_MassWeightedVector_MassWeightedVector(this,&
   & that) result(output)
  implicit none
  
  type(MassWeightedVector), intent(in) :: this
  type(MassWeightedVector), intent(in) :: that
  type(MassWeightedVector)             :: output
  
  output = MassWeightedVector(this%vectors - that%vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MassWeightedVector(this,input)
  implicit none
  
  class(MassWeightedVector), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  select type(this); type is(MassWeightedVector)
    this = MassWeightedVector(RealVector(input))
  end select
end subroutine

function write_MassWeightedVector(this) result(output)
  implicit none
  
  class(MassWeightedVector), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(MassWeightedVector)
    output = str(this%vectors)
  end select
end function

impure elemental function new_MassWeightedVector_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MassWeightedVector)      :: this
  
  this = input
end function
end module
