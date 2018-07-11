! ======================================================================
! A minimal representation of the SymmetryOperator class.
! ======================================================================
module basic_symmetry_submodule
  use utils_module
  implicit none
  
  private
  
  public :: BasicSymmetry
  
  type, extends(Stringsable) :: BasicSymmetry
    integer          :: id
    type(IntMatrix)  :: rotation
    type(RealVector) :: translation
  contains
    procedure, public :: read  => read_BasicSymmetry
    procedure, public :: write => write_BasicSymmetry
  end type
  
  interface BasicSymmetry
    module procedure new_BasicSymmetry
    module procedure new_BasicSymmetry_Strings
    module procedure new_BasicSymmetry_StringArray
  end interface
contains

! Constructor.
function new_BasicSymmetry(id,rotation,translation) result(output)
  implicit none
  
  integer,          intent(in) :: id
  type(IntMatrix),  intent(in) :: rotation
  type(RealVector), intent(in) :: translation
  type(BasicSymmetry)          :: output
  
  output%id          = id
  output%rotation    = rotation
  output%translation = translation
end function

! I/O.
subroutine read_BasicSymmetry(this,input)
  implicit none
  
  class(BasicSymmetry), intent(out) :: this
  type(String),         intent(in)  :: input(:)
  
  integer          :: id
  type(IntMatrix)  :: rotation
  type(RealVector) :: translation
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(BasicSymmetry)
    line = split_line(input(1))
    id = int(line(2))
    rotation = IntMatrix(input(3:5))
    translation = RealVector(input(7))
    
    this = BasicSymmetry(id,rotation,translation)
  class default
    call err()
  end select
end subroutine

function write_BasicSymmetry(this) result(output)
  implicit none
  
  class(BasicSymmetry), intent(in) :: this
  type(String), allocatable        :: output(:)
  
  select type(this); type is(BasicSymmetry)
    output = [ 'Operation '//this%id, &
             & str('Rotation:'),      &
             & str(this%rotation),    &
             & str('Translation'),    &
             & str(this%translation)  ]
  class default
    call err()
  end select
end function

function new_BasicSymmetry_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(BasicSymmetry)      :: this
  
  call this%read(input)
end function

impure elemental function new_BasicSymmetry_StringArray(input) result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(BasicSymmetry)           :: this
  
  this = BasicSymmetry(str(input))
end function
end module
