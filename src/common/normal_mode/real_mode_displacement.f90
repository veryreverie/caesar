! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use mass_weighted_displacement_submodule
  use cartesian_displacement_submodule
  use real_mode_submodule
  use real_single_mode_vector_submodule
  use real_mode_vector_submodule
  implicit none
  
  private
  
  public :: RealModeDisplacement
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(RealModeVector) :: RealModeDisplacement
  contains
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface RealModeDisplacement
    module procedure new_RealModeDisplacement_RealModeVector
    module procedure new_RealModeDisplacement_RealSingleModeVectors
    module procedure new_RealModeDisplacement_MassWeightedDisplacement
    module procedure new_RealModeDisplacement_CartesianDisplacement
    module procedure new_RealModeDisplacement_StringArray
  end interface
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_RealModeDisplacement
  end interface
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_RealModeDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_RealModeDisplacement
    module procedure multiply_RealModeDisplacement_real
  end interface
  
  interface operator(/)
    module procedure divide_RealModeDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_RealModeDisplacement_RealModeDisplacement
  end interface
  
  interface sum
    module procedure sum_RealModeDisplacements
  end interface
  
  interface operator(-)
    module procedure negative_RealModeDisplacement
    module procedure subtract_RealModeDisplacement_RealModeDisplacement
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
! Arithmetic.
! ----------------------------------------------------------------------
impure elemental function multiply_real_RealModeDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                   intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this*that%RealModeVector)
end function

impure elemental function multiply_RealModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%RealModeVector*that)
end function

impure elemental function divide_RealModeDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  real(dp),                   intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%RealModeVector/that)
end function

impure elemental function add_RealModeDisplacement_RealModeDisplacement(this, &
   & that) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%RealModeVector + that%RealModeVector)
end function

function sum_RealModeDisplacements(this) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this(:)
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(sum(this%RealModeVector))
end function

impure elemental function negative_RealModeDisplacement(this) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(-this%RealModeVector)
end function

impure elemental function subtract_RealModeDisplacement_RealModeDisplacement( &
   & this,that) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: this
  type(RealModeDisplacement), intent(in) :: that
  type(RealModeDisplacement)             :: output
  
  output = RealModeDisplacement(this%RealModeVector - that%RealModeVector)
end function

! ----------------------------------------------------------------------
! Conversions between CartesianDisplacement and RealModeDisplacement.
! ----------------------------------------------------------------------
! Returns the displacement in mass-weighted co-ordinates.
function new_MassWeightedDisplacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(MassWeightedDisplacement)          :: output
  
  output = MassWeightedDisplacement( &
     & MassWeightedVector(this,structure,modes,qpoints))
end function

! Returns the displacement in cartesian co-ordinates.
function new_CartesianDisplacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(CartesianVector(this,structure,modes,qpoints))
end function

! Converts a MassWeightedDisplacement to a RealModeDisplacement.
function new_RealModeDisplacement_MassWeightedDisplacement(displacement, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: displacement
  type(StructureData),            intent(in) :: structure
  type(RealMode),                 intent(in) :: modes(:)
  type(QpointData),               intent(in) :: qpoints(:)
  type(RealModeDisplacement)                 :: this
  
  this = RealModeDisplacement( &
     & RealModeVector(displacement,structure,modes,qpoints))
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
