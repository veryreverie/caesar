! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module complex_mode_displacement_submodule
  use utils_module
  
  use complex_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  
  type, extends(Printable) :: ComplexModeDisplacement
    type(ComplexSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: to_String => to_String_ComplexModeDisplacement
  end type
  
  interface size
    module procedure size_ComplexModeDisplacement
  end interface
contains

! Return the number of modes along which the vector has displacements.
function size_ComplexModeDisplacement(input) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: input
  integer                                   :: output
  
  output = size(input%displacements)
end function

! I/O.
function to_String_ComplexModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = str(this%displacements(i))
  enddo
end function
end module
