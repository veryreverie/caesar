! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module mode_displacement_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use single_mode_displacement_module
  use stringable_module
  use printable_module
  implicit none
  
  private
  
  public :: ModeDisplacement
  
  type, extends(Printable) :: ModeDisplacement
    type(SingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure :: str => str_ModeDisplacement
  end type
  
  interface size
    module procedure size_ModeDisplacement
  end interface
contains

! Return the number of modes along which the vector has displacements.
function size_ModeDisplacement(input) result(output)
  implicit none
  
  type(ModeDisplacement), intent(in) :: input
  integer                            :: output
  
  output = size(input%displacements)
end function

! I/O.
function str_ModeDisplacement(this) result(output)
  implicit none
  
  class(ModeDisplacement), intent(in) :: this
  type(String), allocatable           :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = str(this%displacements(i))
  enddo
end function
end module
