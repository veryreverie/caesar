! ======================================================================
! A displacement in cartesian co-ordinates.
! ======================================================================
module cartesian_displacement_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  implicit none
  
  private
  
  public :: CartesianDisplacement
  public :: displace_structure
  public :: operator(*)
  public :: operator(+)
  public :: sum
  
  type, extends(CartesianVector) :: CartesianDisplacement
  contains
    procedure, public :: read  => read_CartesianDisplacement
    procedure, public :: write => write_CartesianDisplacement
  end type
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement_CartesianVector
    module procedure new_CartesianDisplacement
    module procedure new_CartesianDisplacement_StringArray
  end interface
  
  interface operator(*)
    module procedure multiply_real_CartesianDisplacement
    module procedure multiply_CartesianDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_CartesianDisplacement_CartesianDisplacement
  end interface
  
  interface sum
    module procedure sum_CartesianDisplacements
  end interface
contains

! ----------------------------------------------------------------------
! Constructors.
! ----------------------------------------------------------------------
function new_CartesianDisplacement_CartesianVector(displacements) result(this)
  implicit none
  
  type(CartesianVector), intent(in) :: displacements
  type(CartesianDisplacement)       :: this
  
  this%CartesianVector = displacements
end function

function new_CartesianDisplacement(displacements) result(this)
  implicit none
  
  type(RealVector), intent(in) :: displacements(:)
  type(CartesianDisplacement)  :: this
  
  this = CartesianDisplacement(CartesianVector(displacements))
end function

! ----------------------------------------------------------------------
! Construct the structure which is displaced from the input structure.
! ----------------------------------------------------------------------
function displace_structure(structure,displacement) result(output)
  implicit none
  
  type(StructureData),         intent(in) :: structure
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData)                     :: output
  
  integer :: i
  
  if (structure%no_atoms/=size(displacement)) then
    call print_line(CODE_ERROR//': Trying to displace a structure by a &
       &displacement which does not match the number of atoms.')
    call err()
  endif
  
  output = structure
  
  do i=1,size(displacement)
    call output%atoms(i)%set_cartesian_position( &
        &   output%atoms(i)%cartesian_position() &
        & + displacement%vectors(i))
  enddo
end function

! ----------------------------------------------------------------------
! Algebra.
! ----------------------------------------------------------------------
function multiply_real_CartesianDisplacement(this,that) result(output)
  implicit none
  
  real(dp),                    intent(in) :: this
  type(CartesianDisplacement), intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this * that%CartesianVector)
end function

function multiply_CartesianDisplacement_real(this,that) result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  real(dp),                    intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%CartesianVector * that)
end function

function add_CartesianDisplacement_CartesianDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  type(CartesianDisplacement), intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%CartesianVector + that%CartesianVector)
end function

function sum_CartesianDisplacements(input) result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: input(:)
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(sum(input%CartesianVector))
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianDisplacement(this,input)
  implicit none
  
  class(CartesianDisplacement), intent(out) :: this
  type(String),                 intent(in)  :: input(:)
  
  select type(this); type is(CartesianDisplacement)
    this = CartesianDisplacement(CartesianVector(StringArray(input)))
  end select
end subroutine

function write_CartesianDisplacement(this) result(output)
  implicit none
  
  class(CartesianDisplacement), intent(in) :: this
  type(String), allocatable                :: output(:)
  
  select type(this); type is(CartesianDisplacement)
    output = str(this%CartesianVector)
  end select
end function

impure elemental function new_CartesianDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianDisplacement)   :: this
  
  this = input
end function
end module
