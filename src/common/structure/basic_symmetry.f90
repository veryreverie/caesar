! ======================================================================
! A minimal representation of the SymmetryOperator class.
! ======================================================================
module caesar_basic_symmetry_module
  use caesar_utils_module
  
  use caesar_atom_module
  use caesar_spglib_module
  implicit none
  
  private
  
  public :: BasicSymmetry
  
  type, extends(Stringsable) :: BasicSymmetry
    integer                      :: id
    type(IntMatrix)              :: tensor
    type(RealVector)             :: translation
    type(Group)                  :: atom_group
    type(IntVector), allocatable :: rvectors(:)
  contains
    procedure, public :: read  => read_BasicSymmetry
    procedure, public :: write => write_BasicSymmetry
  end type
  
  interface BasicSymmetry
    ! Constructor.
    module function new_BasicSymmetry(id,tensor,translation,atom_group, &
       & rvectors) result(this) 
      integer,          intent(in) :: id
      type(IntMatrix),  intent(in) :: tensor
      type(RealVector), intent(in) :: translation
      type(Group),      intent(in) :: atom_group
      type(IntVector),  intent(in) :: rvectors(:)
      type(BasicSymmetry)          :: this
    end function
  
    ! ----------------------------------------------------------------------
    ! Calculates symmetry data from Spglib symmetries.
    ! ----------------------------------------------------------------------
    module function new_BasicSymmetry_SpglibSymmetries(input,atoms, &
       & symmetry_precision) result(output) 
      type(SpglibSymmetries), intent(in) :: input
      type(AtomData),         intent(in) :: atoms(:)
      real(dp),               intent(in) :: symmetry_precision
      type(BasicSymmetry), allocatable   :: output(:)
    end function
  end interface
  
  interface
    ! I/O.
    module subroutine read_BasicSymmetry(this,input) 
      class(BasicSymmetry), intent(out) :: this
      type(String),         intent(in)  :: input(:)
    end subroutine
  end interface
  
  interface
    module function write_BasicSymmetry(this) result(output) 
      class(BasicSymmetry), intent(in) :: this
      type(String), allocatable        :: output(:)
    end function
  end interface
  
  interface BasicSymmetry
    module function new_BasicSymmetry_Strings(input) result(this) 
      type(String), intent(in) :: input(:)
      type(BasicSymmetry)      :: this
    end function
  
    impure elemental module function new_BasicSymmetry_StringArray(input) &
       & result(this) 
      type(StringArray), intent(in) :: input
      type(BasicSymmetry)           :: this
    end function
  end interface
end module
