! ======================================================================
! A force in cartesian co-ordinates.
! ======================================================================
module cartesian_force_submodule
  use utils_module
  
  use cartesian_vector_submodule
  implicit none
  
  private
  
  public :: CartesianForce
  public :: operator(*)
  public :: operator(+)
  public :: sum
  
  type, extends(CartesianVector) :: CartesianForce
  contains
    procedure, public :: read  => read_CartesianForce
    procedure, public :: write => write_CartesianForce
  end type
  
  interface CartesianForce
    module procedure new_CartesianForce_CartesianVector
    module procedure new_CartesianForce
    module procedure new_CartesianForce_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_real_CartesianForce
    module procedure multiply_CartesianForce_real
  end interface
  
  interface operator(+)
    module procedure add_CartesianForce_CartesianForce
  end interface
  
  interface sum
    module procedure sum_CartesianForces
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_CartesianForce_CartesianVector(forces) result(this)
  implicit none
  
  type(CartesianVector), intent(in) :: forces
  type(CartesianForce)              :: this
  
  this%CartesianVector = forces
end function

function new_CartesianForce(forces) result(this)
  implicit none
  
  type(RealVector), intent(in) :: forces(:)
  type(CartesianForce)         :: this
  
  this = CartesianForce(CartesianVector(forces))
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
function multiply_real_CartesianForce(this,that) result(output)
  implicit none
  
  real(dp),             intent(in) :: this
  type(CartesianForce), intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this * that%CartesianVector)
end function

function multiply_CartesianForce_real(this,that) result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  real(dp),             intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%CartesianVector * that)
end function

function add_CartesianForce_CartesianForce(this,that) &
   & result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: this
  type(CartesianForce), intent(in) :: that
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%CartesianVector + that%CartesianVector)
end function

function sum_CartesianForces(input) result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: input(:)
  type(CartesianForce)             :: output
  
  output = CartesianForce(sum(input%CartesianVector))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianForce(this,input)
  implicit none
  
  class(CartesianForce), intent(out) :: this
  type(String),          intent(in)  :: input(:)
  
  select type(this); type is(CartesianForce)
    this = CartesianForce(CartesianVector(StringArray(input)))
  end select
end subroutine

function write_CartesianForce(this) result(output)
  implicit none
  
  class(CartesianForce), intent(in) :: this
  type(String), allocatable         :: output(:)
  
  select type(this); type is(CartesianForce)
    output = str(this%CartesianVector)
  end select
end function

impure elemental function new_CartesianForce_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianForce)          :: this
  
  this = input
end function
end module
