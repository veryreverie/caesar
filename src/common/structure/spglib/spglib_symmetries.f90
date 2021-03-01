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
    module function new_SpglibSymmetries(spacegroup_number,             &
       & international_symbol,transformation,origin_shift,n_operations, &
       & tensors,translations,n_atoms,pointgroup_symbol) result(this) 
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
    end function
  end interface
  
  interface size
    module function size_SpglibSymmetries(this) result(output) 
      type(SpglibSymmetries), intent(in) :: this
      integer                            :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! I/O.
    ! ----------------------------------------------------------------------
    module function write_SpglibSymmetries(this) result(output) 
      class(SpglibSymmetries), intent(in) :: this
      type(String), allocatable          :: output(:)
    end function
  end interface
end module
