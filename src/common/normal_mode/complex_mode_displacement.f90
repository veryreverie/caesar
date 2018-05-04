! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module complex_mode_displacement_submodule
  use utils_module
  
  use complex_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  
  type, extends(Stringsable) :: ComplexModeDisplacement
    type(ComplexSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
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

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeDisplacement(this,input)
  implicit none
  
  class(ComplexModeDisplacement), intent(out) :: this
  type(String),                   intent(in)  :: input(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexModeDisplacement)
    allocate(this%displacements(size(input)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      this%displacements(i) = input(i)
    enddo
  end select
end subroutine

function write_ComplexModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  integer :: i,ialloc
  
  select type(this); type is(ComplexModeDisplacement)
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = str(this%displacements(i))
    enddo
  end select
end function
end module
