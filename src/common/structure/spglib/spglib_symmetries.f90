! ======================================================================
! A class to hold the output of Spglib.
! ======================================================================
module spglib_symmetries_module
  use utils_module
  implicit none
  
  private
  
  public :: SpglibSymmetries
  public :: size
  
  type, extends(NoDefaultConstructor) :: SpglibSymmetries
    integer                       :: spacegroup_number
    integer                       :: hall_number
    type(String)                  :: international_symbol
    type(String)                  :: hall_symbol
    type(String)                  :: choice
    type(RealMatrix)              :: transformation
    type(RealVector)              :: origin_shift
    integer                       :: n_operations
    type(IntMatrix),  allocatable :: tensors(:)
    type(RealVector), allocatable :: translations(:)
    integer                       :: n_atoms
    type(String)                  :: pointgroup_symbol
  end type
  
  interface SpglibSymmetries
    module procedure new_SpglibSymmetries
  end interface
  
  interface size
    module procedure size_SpglibSymmetries
  end interface
contains

function new_SpglibSymmetries(spacegroup_number,hall_number,              &
   & international_symbol,hall_symbol,choice,transformation,origin_shift, &
   & n_operations,tensors,translations,n_atoms,pointgroup_symbol) result(this)
  implicit none
  
  integer,          intent(in) :: spacegroup_number
  integer,          intent(in) :: hall_number
  type(String),     intent(in) :: international_symbol
  type(String),     intent(in) :: hall_symbol
  type(String),     intent(in) :: choice
  type(RealMatrix), intent(in) :: transformation
  type(RealVector), intent(in) :: origin_shift
  integer,          intent(in) :: n_operations
  type(IntMatrix),  intent(in) :: tensors(:)
  type(RealVector), intent(in) :: translations(:)
  integer,          intent(in) :: n_atoms
  type(String),     intent(in) :: pointgroup_symbol
  type(SpglibSymmetries)       :: this
  
  if (size(tensors)/=n_operations) then
    call print_line(ERROR//': n_operations does not match the number of &
       &tensors.')
    call err()
  elseif (size(translations)/=n_operations) then
    call print_line(ERROR//': n_operations does not match the number of &
       &translations.')
    call err()
  endif
  
  this%spacegroup_number = spacegroup_number
  this%hall_number = hall_number
  this%choice = choice
  this%transformation = transformation
  this%origin_shift = origin_shift
  this%n_operations = n_operations
  this%tensors = tensors
  this%translations = translations
  this%n_atoms = n_atoms
  this%pointgroup_symbol = pointgroup_symbol
end function

function size_SpglibSymmetries(this) result(output)
  implicit none
  
  type(SpglibSymmetries), intent(in) :: this
  integer                            :: output
  
  output = size(this%tensors)
end function
end module
