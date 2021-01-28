! ======================================================================
! A force in mass-weighted cartesian co-ordinates.
! ======================================================================
module caesar_mass_weighted_force_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_cartesian_force_module
  implicit none
  
  private
  
  public :: MassWeightedForce
  public :: size
  public :: CartesianForce
  public :: operator(*)
  public :: operator(/)
  public :: operator(+)
  public :: operator(-)
  public :: sum
  
  type, extends(Stringsable) :: MassWeightedForce
    type(RealVector), allocatable :: vectors(:)
  contains
    procedure, public :: read  => read_MassWeightedForce
    procedure, public :: write => write_MassWeightedForce
  end type
  
  interface MassWeightedForce
    module procedure new_MassWeightedForce
    module procedure new_MassWeightedForce_zero
    module procedure new_MassWeightedForce_CartesianForce
    module procedure new_MassWeightedForce_Strings
    module procedure new_MassWeightedForce_StringArray
  end interface
  
  interface size
    module procedure size_MassWeightedForce
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
! Constructor and size() function.
! ----------------------------------------------------------------------
function new_MassWeightedForce(forces) result(this)
  implicit none
  
  type(RealVector), intent(in) :: forces(:)
  type(MassWeightedForce)      :: this
  
  this%vectors = forces
end function

function size_MassWeightedForce(this) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  integer                             :: output
  
  output = size(this%vectors)
end function

! ----------------------------------------------------------------------
! Construct a zero force.
! ----------------------------------------------------------------------
impure elemental function new_MassWeightedForce_zero(structure) &
   & result(this)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(MassWeightedForce)         :: this
  
  integer :: i
  
  this%vectors = [(dblevec(zeroes(3)), i=1, structure%no_atoms)]
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
  
  output = MassWeightedForce(input%vectors/sqrt(structure%atoms%mass()))
end function

function new_CartesianForce_MassWeightedForce(input,structure) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: input
  type(StructureData),     intent(in) :: structure
  type(CartesianForce)                :: output
  
  output = CartesianForce(input%vectors*sqrt(structure%atoms%mass()))
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
  
  output = MassWeightedForce(this * that%vectors)
end function

impure elemental function multiply_MassWeightedForce_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  real(dp),                intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%vectors * that)
end function

impure elemental function divide_MassWeightedForce_real(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  real(dp),                intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%vectors / that)
end function

impure elemental function add_MassWeightedForce_MassWeightedForce(this,that) &
   & result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce), intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%vectors + that%vectors)
end function

function sum_MassWeightedForces(input) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: input(:)
  type(MassWeightedForce)             :: output
  
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

impure elemental function negative_MassWeightedForce(this) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(-this%vectors)
end function

impure elemental function subtract_MassWeightedForce_MassWeightedForce(this, &
   & that) result(output)
  implicit none
  
  type(MassWeightedForce), intent(in) :: this
  type(MassWeightedForce), intent(in) :: that
  type(MassWeightedForce)             :: output
  
  output = MassWeightedForce(this%vectors - that%vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_MassWeightedForce(this,input)
  implicit none
  
  class(MassWeightedForce), intent(out) :: this
  type(String),             intent(in)  :: input(:)
  
  select type(this); type is(MassWeightedForce)
    this = MassWeightedForce(RealVector(input))
  class default
    call err()
  end select
end subroutine

function write_MassWeightedForce(this) result(output)
  implicit none
  
  class(MassWeightedForce), intent(in) :: this
  type(String), allocatable            :: output(:)
  
  select type(this); type is(MassWeightedForce)
    output = str(this%vectors)
  class default
    call err()
  end select
end function

function new_MassWeightedForce_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(MassWeightedForce)  :: this
  
  call this%read(input)
end function

impure elemental function new_MassWeightedForce_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(MassWeightedForce)       :: this
  
  this = MassWeightedForce(str(input))
end function
end module
