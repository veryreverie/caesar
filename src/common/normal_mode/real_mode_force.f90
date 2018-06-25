! ======================================================================
! A force in real mode co-ordinates.
! ======================================================================
module real_mode_force_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_force_submodule
  use mass_weighted_force_submodule
  use real_mode_submodule
  use real_single_mode_vector_submodule
  use real_mode_vector_submodule
  implicit none
  
  private
  
  public :: RealModeForce
  public :: MassWeightedForce
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(RealModeVector) :: RealModeForce
  contains
    procedure, public :: read  => read_RealModeForce
    procedure, public :: write => write_RealModeForce
  end type
  
  interface RealModeForce
    module procedure new_RealModeForce_RealModeVector
    module procedure new_RealModeForce_RealSingleModeVectors
    module procedure new_RealModeForce_RealModes
    module procedure new_RealModeForce_MassWeightedForce
    module procedure new_RealModeForce_CartesianForce
    module procedure new_RealModeForce_StringArray
  end interface
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce_RealModeForce
  end interface
  
  interface CartesianForce
    module procedure new_CartesianForce_RealModeForce
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealModeForce
    module procedure multiply_RealModeForce_real
  end interface
  
  interface operator(/)
    module procedure divide_RealModeForce_real
  end interface
  
  interface operator(+)
    module procedure add_RealModeForce_RealModeForce
  end interface
  
  interface sum
    module procedure sum_RealModeForces
  end interface
  
  interface operator(-)
    module procedure negative_RealModeForce
    module procedure subtract_RealModeForce_RealModeForce
  end interface
contains

! Constructors.
function new_RealModeForce_RealModeVector(force) result(this)
  implicit none
  
  type(RealModeVector), intent(in) :: force
  type(RealModeForce)              :: this
  
  this%RealModeVector = force
end function

function new_RealModeForce_RealSingleModeVectors(forces) &
   & result(this)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: forces(:)
  type(RealModeForce)                    :: this
  
  this = RealModeForce(RealModeVector(forces))
end function

function new_RealModeForce_RealModes(modes,forces) result(this)
  implicit none
  
  type(RealMode), intent(in) :: modes(:)
  real(dp),       intent(in) :: forces(:)
  type(RealModeForce)        :: this
  
  this = RealModeForce(RealModeVector(modes,forces))
end function

! ----------------------------------------------------------------------
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealModeForce(this,that) &
   & result(output)
  implicit none
  
  real(dp),            intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this*that%RealModeVector)
end function

impure elemental function multiply_RealModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  real(dp),            intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%RealModeVector*that)
end function

impure elemental function divide_RealModeForce_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  real(dp),            intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%RealModeVector/that)
end function

impure elemental function add_RealModeForce_RealModeForce(this, &
   & that) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%RealModeVector + that%RealModeVector)
end function

function sum_RealModeForces(this) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this(:)
  type(RealModeForce)             :: output
  
  output = RealModeForce(sum(this%RealModeVector))
end function

impure elemental function negative_RealModeForce(this) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce)             :: output
  
  output = RealModeForce(-this%RealModeVector)
end function

impure elemental function subtract_RealModeForce_RealModeForce( &
   & this,that) result(output)
  implicit none
  
  type(RealModeForce), intent(in) :: this
  type(RealModeForce), intent(in) :: that
  type(RealModeForce)             :: output
  
  output = RealModeForce(this%RealModeVector - that%RealModeVector)
end function

! ----------------------------------------------------------------------
! Conversions between CartesianForce and RealModeForce.
! ----------------------------------------------------------------------
! Returns the force in mass-weighted co-ordinates.
function new_MassWeightedForce_RealModeForce(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(MassWeightedForce)          :: output
  
  output = MassWeightedForce(MassWeightedVector(this,structure,modes,qpoints))
end function

! Returns the force in cartesian co-ordinates.
function new_CartesianForce_RealModeForce(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(CartesianForce)             :: output
  
  output = CartesianForce(CartesianVector(this,structure,modes,qpoints))
end function

! Converts a MassWeightedForce to a RealModeForce.
function new_RealModeForce_MassWeightedForce(force, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(MassWeightedForce), intent(in) :: force
  type(StructureData),     intent(in) :: structure
  type(RealMode),          intent(in) :: modes(:)
  type(QpointData),        intent(in) :: qpoints(:)
  type(RealModeForce)                 :: this
  
  this = RealModeForce(RealModeVector(force,structure,modes,qpoints))
end function

! Converts a CartesianForce to a RealModeForce.
function new_RealModeForce_CartesianForce(force, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(CartesianForce), intent(in) :: force
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(RealModeForce)              :: this
  
  this = RealModeForce(RealModeVector(force,structure,modes,qpoints))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeForce(this,input)
  implicit none
  
  class(RealModeForce), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  select type(this); type is(RealModeForce)
    this = RealModeForce(RealModeVector(StringArray(input)))
  end select
end subroutine

function write_RealModeForce(this) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(RealModeForce)
    output = str(this%RealModeVector)
  end select
end function

impure elemental function new_RealModeForce_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeForce)    :: this
  
  this = input
end function
end module
