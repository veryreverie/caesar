! ======================================================================
! A force in mass-weighted cartesian co-ordinates.
! ======================================================================
module mass_weighted_force_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  use cartesian_force_submodule
  use mass_weighted_vector_submodule
  implicit none
  
  private
  
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  
  type, extends(MassWeightedVector) :: MassWeightedForce
  contains
    procedure, public :: read  => read_MassWeightedForce
    procedure, public :: write => write_MassWeightedForce
  end type
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce_MassWeightedVector
    module procedure new_MassWeightedForce
    module procedure new_MassWeightedForce_CartesianForce
    module procedure new_MassWeightedForce_StringArray
  end interface
  
  interface CartesianForce
    module procedure new_CartesianForce_MassWeightedForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_MassWeightedForce
    module procedure multiply_MassWeightedForce_real
  end interface
  
  interface operator(/)
    module procedure divide_MassWeightedForce_real
  end interface
  
  interface operator(+)
    module procedure add_MassWeightedForce_MassWeightedForce
  end interface
  
  interface sum
    module procedure sum_MassWeightedForces
  end interface
  
  interface operator(-)
    module procedure negative_MassWeightedForce
    module procedure subtract_MassWeightedForce_MassWeightedForce
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_MassWeightedForce_MassWeightedVector(forces) result(this)
  implicit none
  
  type(MassWeightedVector), intent(in) :: forces
  type(MassWeightedForce)              :: this
  
  this%MassWeightedVector = forces
end function

function new_MassWeightedForce(forces) result(this)
  implicit none
  
  type(RealVector), intent(in) :: forces(:)
  type(MassWeightedForce)      :: this
  
  this = MassWeightedForce(MassWeightedVector(forces))
end function

! ----------------------------------------------------------------------
! Conversion to and from non-mass-weighted co-ordinates.
! ----------------------------------------------------------------------
function new_MassWeightedForce_CartesianForce(input,structure) &
   & result(output)
  implicit none
  
  type(CartesianForce), intent(in) :: input
  type(StructureData),  intent(in) :: structure
  type(MassWeightedForce)          :: output
  
  output = MassWeightedForce(MassWeightedVector(input,structure))
end function

function new_CartesianForce_MassWeightedForce(input,structure) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: input
  type(StructureData),     intent(in) :: structure
  type(CartesianForce)                :: output
  
  output = CartesianForce(CartesianVector(input,structure))
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
impure elemental function multiply_real_MassWeightedForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),                intent(in) :: this
  type(MassWeightedForce), intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this * that%MassWeightedVector)
end function

impure elemental function multiply_MassWeightedForce_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  real(dp),                intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%MassWeightedVector * that)
end function

impure elemental function divide_MassWeightedForce_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  real(dp),                intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%MassWeightedVector / that)
end function

impure elemental function add_MassWeightedForce_MassWeightedForce(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce), intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%MassWeightedVector + that%MassWeightedVector)
end function

function sum_MassWeightedForces(input) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: input(:)
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(sum(input%MassWeightedVector))
end function

impure elemental function negative_MassWeightedForce(this) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(-this%MassWeightedVector)
end function

impure elemental function subtract_MassWeightedForce_MassWeightedForce(this, &
   & that) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce), intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%MassWeightedVector - that%MassWeightedVector)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MassWeightedForce(this,input)
  implicit none
  
  class(MassWeightedForce), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  select type(this); type is(MassWeightedForce)
    this = MassWeightedForce(MassWeightedVector(StringArray(input)))
  end select
end subroutine

function write_MassWeightedForce(this) result(output)
  implicit none
  
  class(MassWeightedForce), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(MassWeightedForce)
    output = str(this%MassWeightedVector)
  end select
end function

impure elemental function new_MassWeightedForce_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MassWeightedForce)       :: this
  
  this = input
end function
end module
