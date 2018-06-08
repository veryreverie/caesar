! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_displacement_submodule
  use real_mode_submodule
  use real_single_mode_vector_submodule
  use real_mode_vector_submodule
  implicit none
  
  private
  
  public :: RealModeDisplacement
  
  type, extends(RealModeVector) :: RealModeDisplacement
  contains
    procedure, public :: cartesian_displacement => &
       & cartesian_displacement_RealModeDisplacement
    
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface RealModeDisplacement
    module procedure new_RealModeDisplacement_RealModeVector
    module procedure new_RealModeDisplacement_RealSingleModeVectors
    module procedure new_RealModeDisplacement_CartesianDisplacement
    module procedure new_RealModeDisplacement_StringArray
  end interface
contains

! Constructors.
function new_RealModeDisplacement_RealModeVector(displacement) result(this)
  implicit none
  
  type(RealModeVector), intent(in) :: displacement
  type(RealModeDisplacement)       :: this
  
  this%RealModeVector = displacement
end function

function new_RealModeDisplacement_RealSingleModeVectors(displacements) &
   & result(this)
  implicit none
  
  type(RealSingleModeVector), intent(in) :: displacements(:)
  type(RealModeDisplacement)             :: this
  
  this = RealModeDisplacement(RealModeVector(displacements))
end function

! ----------------------------------------------------------------------
! Conversions between CartesianDisplacement and RealModeDisplacement.
! ----------------------------------------------------------------------
! Returns the displacement in cartesian co-ordinates.
function cartesian_displacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement( &
     & this%cartesian_vector(structure,modes,qpoints))
end function

! Converts a CartesianDisplacement to a RealModeDisplacement.
function new_RealModeDisplacement_CartesianDisplacement(displacement, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(RealModeDisplacement)              :: this
  
  this = RealModeDisplacement( &
     & RealModeVector(displacement,structure,modes,qpoints))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeDisplacement(this,input)
  implicit none
  
  class(RealModeDisplacement), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  select type(this); type is(RealModeDisplacement)
    this = RealModeDisplacement(RealModeVector(StringArray(input)))
  end select
end subroutine

function write_RealModeDisplacement(this) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(RealModeDisplacement)
    output = str(this%RealModeVector)
  end select
end function

impure elemental function new_RealModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeDisplacement)    :: this
  
  this = input
end function
end module
