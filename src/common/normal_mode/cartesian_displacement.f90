! ======================================================================
! A displacement in cartesian co-ordinates.
! ======================================================================
module cartesian_displacement_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: CartesianDisplacement
  public :: size
  public :: displace_structure
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: sum
  public :: operator(-)
  
  type, extends(Stringsable) :: CartesianDisplacement
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_CartesianDisplacement
    procedure, public :: write => write_CartesianDisplacement
  end type
  
  interface CartesianDisplacement
    module procedure new_CartesianDisplacement
    module procedure new_CartesianDisplacement_Strings
    module procedure new_CartesianDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_CartesianDisplacement
  end interface
  
  interface displace_structure
    module procedure displace_structure_CartesianDisplacement
  end interface
  
  interface operator(*)
    module procedure multiply_real_CartesianDisplacement
    module procedure multiply_CartesianDisplacement_real
  end interface
  
  interface operator(/)
    module procedure divide_CartesianDisplacement_real
  end interface
  
  interface operator(+)
    module procedure add_CartesianDisplacement_CartesianDisplacement
  end interface
  
  interface sum
    module procedure sum_CartesianDisplacements
  end interface
  
  interface operator(-)
    module procedure negative_CartesianDisplacement
    module procedure subtract_CartesianDisplacement_CartesianDisplacement
  end interface
contains

! ----------------------------------------------------------------------
! Constructor and size() function.
! ----------------------------------------------------------------------
function new_CartesianDisplacement(displacements) result(this)
  implicit none
  
  type(RealVector), intent(in) :: displacements(:)
  type(CartesianDisplacement)  :: this
  
  this%vectors = displacements
end function

function size_CartesianDisplacement(this) result(output)
  implicit none
  
  class(CartesianDisplacement), intent(in) :: this
  integer                                  :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Construct the structure which is displaced from the input structure.
! ----------------------------------------------------------------------
function displace_structure_CartesianDisplacement(structure,displacement) &
   & result(output)
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
impure elemental function multiply_real_CartesianDisplacement(this,that) &
   & result(output)
  implicit none
  
  real(dp),                    intent(in) :: this
  type(CartesianDisplacement), intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this * that%vectors)
end function

impure elemental function multiply_CartesianDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  real(dp),                    intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%vectors * that)
end function

impure elemental function divide_CartesianDisplacement_real(this,that) &
   & result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  real(dp),                    intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%vectors / that)
end function

impure elemental function add_CartesianDisplacement_CartesianDisplacement( &
   & this,that) result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  type(CartesianDisplacement), intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%vectors + that%vectors)
end function

function sum_CartesianDisplacements(input) result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: input(:)
  type(CartesianDisplacement)             :: output
  
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

impure elemental function negative_CartesianDisplacement(this) result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(-this%vectors)
end function

impure elemental function                                            &
   & subtract_CartesianDisplacement_CartesianDisplacement(this,that) &
   & result(output)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: this
  type(CartesianDisplacement), intent(in) :: that
  type(CartesianDisplacement)             :: output
  
  output = CartesianDisplacement(this%vectors - that%vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_CartesianDisplacement(this,input)
  implicit none
  
  class(CartesianDisplacement), intent(out) :: this
  type(String),                 intent(in)  :: input(:)
  
  select type(this); type is(CartesianDisplacement)
    this = CartesianDisplacement(RealVector(input))
  class default
    call err()
  end select
end subroutine

function write_CartesianDisplacement(this) result(output)
  implicit none
  
  class(CartesianDisplacement), intent(in) :: this
  type(String), allocatable                :: output(:)
  
  select type(this); type is(CartesianDisplacement)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_CartesianDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)    :: input(:)
  type(CartesianDisplacement) :: this
  
  call this%read(input)
end function

impure elemental function new_CartesianDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(CartesianDisplacement)   :: this
  
  this = CartesianDisplacement(str(input))
end function
end module
