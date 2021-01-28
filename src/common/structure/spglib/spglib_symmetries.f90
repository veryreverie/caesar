! ======================================================================
! A class to hold the output of Spglib.
! ======================================================================
module caesar_spglib_symmetries_module
  use caesar_utils_module
  implicit none
  
  private
  
  public :: SpglibSymmetries
  public :: size
  
  type, extends(Stringswriteable) :: SpglibSymmetries
    integer                       :: spacegroup_number
    type(String)                  :: international_symbol
    type(RealMatrix)              :: transformation
    type(RealVector)              :: origin_shift
    integer                       :: n_operations
    type(IntMatrix),  allocatable :: tensors(:)
    type(RealVector), allocatable :: translations(:)
    integer                       :: n_atoms
    type(String)                  :: pointgroup_symbol
  contains
    procedure, public :: write => write_SpglibSymmetries
  end type
  
  interface SpglibSymmetries
    module procedure new_SpglibSymmetries
  end interface
  
  interface size
    module procedure size_SpglibSymmetries
  end interface
contains

function new_SpglibSymmetries(spacegroup_number,international_symbol,       &
   & transformation,origin_shift,n_operations,tensors,translations,n_atoms, &
   & pointgroup_symbol) result(this)
  implicit none
  
  integer,          intent(in) :: spacegroup_number
  type(String),     intent(in) :: international_symbol
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
  this%international_symbol = international_symbol
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

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
function write_SpglibSymmetries(this) result(output)
  implicit none
  
  class(SpglibSymmetries), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  integer :: i
  
  select type(this); type is(SpglibSymmetries)
    output = [ 'Spacegroup Number    : '//this%spacegroup_number,    &
             & 'International Symbol : '//this%international_symbol, &
             & 'Pointgroup Symbol    : '//this%pointgroup_symbol,    &
             & str('Transformation       : '),                       &
             & str(this%transformation),                             &
             & 'Origin Shift         : '//this%origin_shift,         &
             & 'No. Atoms            : '//this%n_atoms,              &
             & 'No. Operations       : '//this%n_operations,         &
             & str('Operations           : ')                        ]
    do i=1,this%n_operations
      output = [ output,                   &
               & str(this%tensors(i)),     &
               & str(this%translations(i)) ]
    enddo
  class default
    call err()
  end select
end function
end module
