! ======================================================================
! A displacement in mass-weighted cartesian co-ordinates.
! ======================================================================
module caesar_mass_weighted_displacement_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_displacement_module
  implicit none
  
  private
  
  public :: MassWeightedDisplacement
  public :: size
  public :: CartesianDisplacement
  public :: displace_structure
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: MassWeightedDisplacement
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_MassWeightedDisplacement
    procedure, public :: write => write_MassWeightedDisplacement
  end type
  
  interface MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement
    module procedure new_MassWeightedDisplacement_zero
    module procedure new_MassWeightedDisplacement_CartesianDisplacement
    module procedure new_MassWeightedDisplacement_Strings
    module procedure new_MassWeightedDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_MassWeightedDisplacement
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
! Constructor and size() function.
! ----------------------------------------------------------------------
function new_MassWeightedDisplacement(displacements) result(this)
  implicit none
  
  type(RealVector), intent(in)   :: displacements(:)
  type(MassWeightedDisplacement) :: this
  
  this%vectors = displacements
end function

function size_MassWeightedDisplacement(this) result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  integer                                    :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Construct a zero displacement.
! ----------------------------------------------------------------------
impure elemental function new_MassWeightedDisplacement_zero(structure) &
   & result(this)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(MassWeightedDisplacement)  :: this
  
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
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
  
  output = MassWeightedDisplacement(input%vectors*sqrt(structure%atoms%mass()))
end function

function new_CartesianDisplacement_MassWeightedDisplacement(input,structure) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: input
  type(StructureData),            intent(in) :: structure
  type(CartesianDisplacement)                :: output
  
  output = CartesianDisplacement(input%vectors/sqrt(structure%atoms%mass()))
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
  
  output = MassWeightedDisplacement(this * that%vectors)
end function

impure elemental function multiply_MassWeightedDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  real(dp),                       intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%vectors * that)
end function

impure elemental function divide_MassWeightedDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  real(dp),                       intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%vectors / that)
end function

impure elemental function                                             &
   & add_MassWeightedDisplacement_MassWeightedDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement), intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%vectors + that%vectors)
end function

function sum_MassWeightedDisplacements(input) result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: input(:)
  type(MassWeightedDisplacement)             :: output
  
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

impure elemental function negative_MassWeightedDisplacement(this) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(-this%vectors)
end function

impure elemental function                                                  &
   & subtract_MassWeightedDisplacement_MassWeightedDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedDisplacement), intent(in) :: this
  type(MassWeightedDisplacement), intent(in) :: that
  type(MassWeightedDisplacement)             :: output
  
  output = MassWeightedDisplacement(this%vectors - that%vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MassWeightedDisplacement(this,input)
  implicit none
  
  class(MassWeightedDisplacement), intent(out) :: this
  type(String),                    intent(in)  :: input(:)
  
  select type(this); type is(MassWeightedDisplacement)
    this = MassWeightedDisplacement(RealVector(input))
  class default
    call err()
  end select
end subroutine

function write_MassWeightedDisplacement(this) result(output)
  implicit none
  
  class(MassWeightedDisplacement), intent(in) :: this
  type(String), allocatable                   :: output(:)
  
  select type(this); type is(MassWeightedDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_MassWeightedDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)       :: input(:)
  type(MassWeightedDisplacement) :: this
  
  call this%read(input)
end function

impure elemental function new_MassWeightedDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in)  :: input
  type(MassWeightedDisplacement) :: this
  
  this = MassWeightedDisplacement(str(input))
end function
end module
