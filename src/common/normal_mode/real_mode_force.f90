! ======================================================================
! A force in real mode co-ordinates.
! ======================================================================
module real_mode_force_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_force_submodule
  use real_mode_submodule
  use real_single_mode_vector_submodule
  use real_mode_vector_submodule
  implicit none
  
  private
  
  public :: RealModeForce
  
  type, extends(RealModeVector) :: RealModeForce
  contains
    procedure, public :: cartesian_force => cartesian_force_RealModeForce
    
    procedure, public :: read  => read_RealModeForce
    procedure, public :: write => write_RealModeForce
  end type
  
  interface RealModeForce
    module procedure new_RealModeForce_RealModeVector
    module procedure new_RealModeForce_RealSingleModeVectors
    module procedure new_RealModeForce_CartesianForce
    module procedure new_RealModeForce_StringArray
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

! ----------------------------------------------------------------------
! Conversions between CartesianForce and RealModeForce.
! ----------------------------------------------------------------------
! Returns the force in cartesian co-ordinates.
function cartesian_force_RealModeForce(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeForce), intent(in) :: this
  type(StructureData),  intent(in) :: structure
  type(RealMode),       intent(in) :: modes(:)
  type(QpointData),     intent(in) :: qpoints(:)
  type(CartesianForce)             :: output
  
  output = CartesianForce(this%cartesian_vector(structure,modes,qpoints))
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
