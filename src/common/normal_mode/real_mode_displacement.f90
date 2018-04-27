! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_submodule
  use utils_module
  
  use real_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: RealModeDisplacement
  
  type, extends(Printable) :: RealModeDisplacement
    type(RealSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: to_String => to_String_RealModeDisplacement
  end type
  
  interface size
    module procedure size_RealModeDisplacement
  end interface
  
  interface RealModeDisplacement
    module procedure new_RealModeDisplacement_Strings
  end interface
contains

! Return the number of modes along which the vector has displacements.
function size_RealModeDisplacement(input) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: input
  integer                                :: output
  
  output = size(input%displacements)
end function

! I/O.
function to_String_RealModeDisplacement(this) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = str(this%displacements(i))
  enddo
end function

function new_RealModeDisplacement_Strings(input) result(this)
  implicit none
  
  type(String), intent(in)   :: input(:)
  type(RealModeDisplacement) :: this
  
  integer :: i,ialloc
  
  allocate(this%displacements(size(input)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    this%displacements(i) = RealSingleModeDisplacement(input(i))
  enddo
end function
end module
