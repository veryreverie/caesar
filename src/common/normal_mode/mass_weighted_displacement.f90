! ======================================================================
! A displacement in mass-weighted cartesian co-ordinates.
! ======================================================================
module mass_weighted_displacement_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  use cartesian_displacement_submodule
  use mass_weighted_vector_submodule
  implicit none
  
  private
  
  public :: MassWeightedDisplacement
  public :: CartesianDisplacement
  public :: displace_structure
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(MassWeightedVector) :: MassWeightedDisplacement
  contains
    procedure, public :: read  => read_MassWeightedDisplacement
    procedure, public :: write => write_MassWeightedDisplacement
  end type
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_MassWeightedVector
    module procedure new_MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_CartesianDisplacement
    module procedure new_MassWeightedDisplacement_StringArray
  end interface
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_MassWeightedDisplacement
  end interface
  
  interface displace_structure
    module procedure displace_structure_MassWeightedDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_MassWeightedDisplacement
    module procedure multiply_MassWeightedDisplacement_real
  end interface
  
  interface operator(/)
    module procedure divide_MassWeightedDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_MassWeightedDisplacement_MassWeightedDisplacement
  end interface
  
  interface sum
    module procedure sum_MassWeightedDisplacements
  end interface
  
  interface operator(-)
    module procedure displace_structure_MassWeightedDisplacement
    module procedure negative_MassWeightedDisplacement
    module procedure subtract_MassWeightedDisplacement_MassWeightedDisplacement
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_MassWeightedDisplacement_MassWeightedVector(displacements) &
   & result(this)
  implicit none
  
  type(MassWeightedVector), intent(in) :: displacements
  type(MassWeightedDisplacement)       :: this
  
  this%MassWeightedVector = displacements
end function

function new_MassWeightedDisplacement(displacements) result(this)
  implicit none
  
  type(RealVector), intent(in)   :: displacements(:)
  type(MassWeightedDisplacement) :: this
  
  this = MassWeightedDisplacement(MassWeightedVector(displacements))
end function

! ----------------------------------------------------------------------
! Conversion to and from non-mass-weighted co-ordinates.
! ----------------------------------------------------------------------
function new_MassWeightedDisplacement_CartesianDisplacement(input,structure) &
   & result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: input
  type(StructureData),         intent(in) :: structure
  type(MassWeightedDisplacement)          :: output
  
  output = MassWeightedDisplacement(MassWeightedVector(input,structure))
end function

function new_CartesianDisplacement_MassWeightedDisplacement(input,structure) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: input
  type(StructureData),            intent(in) :: structure
  type(CartesianDisplacement)                :: output
  
  type(CartesianVector) :: temp
  
  output = CartesianDisplacement(CartesianVector(input,structure))
end function

! ----------------------------------------------------------------------
! Construct the structure which is displaced from the input structure.
! ----------------------------------------------------------------------
function displace_structure_MassWeightedDisplacement(structure,displacement) &
   & result(output)
  implicit none
  
  type(StructureData),            intent(in) :: structure
  type(MassWeightedDisplacement), intent(in) :: displacement
  type(StructureData)                        :: output
  
  type(CartesianDisplacement) :: cartesian_displacement
  
  cartesian_displacement = CartesianDisplacement(displacement,structure)
  output = displace_structure(structure,cartesian_displacement)
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
impure elemental function multiply_real_MassWeightedDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                       intent(in) :: this
  type(MassWeightedDisplacement), intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this * that%MassWeightedVector)
end function

impure elemental function multiply_MassWeightedDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  real(dp),                       intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%MassWeightedVector * that)
end function

impure elemental function divide_MassWeightedDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  real(dp),                       intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%MassWeightedVector / that)
end function

impure elemental function                                             &
   & add_MassWeightedDisplacement_MassWeightedDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement), intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement( this%MassWeightedVector &
                                 & + that%MassWeightedVector)
end function

function sum_MassWeightedDisplacements(input) result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: input(:)
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(sum(input%MassWeightedVector))
end function

impure elemental function negative_MassWeightedDisplacement(this) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(-this%MassWeightedVector)
end function

impure elemental function                                                  &
   & subtract_MassWeightedDisplacement_MassWeightedDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement), intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement( this%MassWeightedVector &
                                 & - that%MassWeightedVector)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MassWeightedDisplacement(this,input)
  implicit none
  
  class(MassWeightedDisplacement), intent(out) :: this
  type(String),                    intent(in)  :: input(:)
  
  select type(this); type is(MassWeightedDisplacement)
    this = MassWeightedDisplacement(MassWeightedVector(StringArray(input)))
  end select
end subroutine

function write_MassWeightedDisplacement(this) result(output)
  implicit none
  
  class(MassWeightedDisplacement), intent(in) :: this
  type(String), allocatable                   :: output(:)
  
  select type(this); type is(MassWeightedDisplacement)
    output = str(this%MassWeightedVector)
  end select
end function

impure elemental function new_MassWeightedDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)  :: input
  type(MassWeightedDisplacement) :: this
  
  this = input
end function
end module
