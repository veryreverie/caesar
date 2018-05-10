! ======================================================================
! A displacement along a single complex mode.
! ======================================================================
module complex_single_mode_displacement_submodule
  use utils_module
  implicit none
  
  private
  
  public :: ComplexSingleModeDisplacement
  
  ! The displacement along a single complex mode.
  type, extends(Stringable) :: ComplexSingleModeDisplacement
    ! The id of the mode.
    integer :: id
    
    ! The displacement along the mode.
    complex(dp) :: displacement
  contains
    procedure, public :: read  => read_ComplexSingleModeDisplacement
    procedure, public :: write => write_ComplexSingleModeDisplacement
  end type
  
  interface ComplexSingleModeDisplacement
    module procedure new_ComplexSingleModeDisplacement
  end interface
contains

! Constructor.
function new_ComplexSingleModeDisplacement(id,displacement) result(this)
  implicit none
  
  integer,     intent(in)             :: id
  complex(dp), intent(in)             :: displacement
  type(ComplexSingleModeDisplacement) :: this
  
  this%id           = id
  this%displacement = displacement
end function

! I/O.
subroutine read_ComplexSingleModeDisplacement(this,input)
  implicit none
  
  class(ComplexSingleModeDisplacement), intent(out) :: this
  type(String),                         intent(in)  :: input
  
  type(String), allocatable :: split_string(:)
  integer                   :: id
  complex(dp)               :: displacement
  
  select type(this); type is(ComplexSingleModeDisplacement)
    split_string = split_line(input)
    if (size(split_string)/=3) then
      call print_line(ERROR//': unable to parse complex single mode &
         &displacement from string: '//input)
      call err()
    endif
    
    ! If e.g. id=3 and power=2.1+1.2i then split_string = ["u3","=","2.1+1.2i"]
    ! The 'u' needs stripping off the first element to give the id.
    id = int(slice(split_string(1),2,len(split_string(1))))
    displacement = cmplx(split_string(3))
    
    this = ComplexSingleModeDisplacement(id,displacement)
  end select
end subroutine

function write_ComplexSingleModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexSingleModeDisplacement), intent(in) :: this
  type(String)                                     :: output
  
  select type(this); type is(ComplexSingleModeDisplacement)
    output = 'u'//this%id//' = '//this%displacement
  end select
end function
end module
