! ======================================================================
! A displacement in complex mode co-ordinates.
! ======================================================================
module complex_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use complex_mode_submodule
  use complex_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: ComplexModeDisplacement
  public :: size
  
  type, extends(Stringsable) :: ComplexModeDisplacement
    type(ComplexSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: qpoints => qpoints_ComplexModeDisplacement
    
    procedure, public :: read  => read_ComplexModeDisplacement
    procedure, public :: write => write_ComplexModeDisplacement
  end type
  
  interface ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement
    module procedure new_ComplexModeDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_ComplexModeDisplacement
  end interface
contains

! Constructor.
function new_ComplexModeDisplacement(displacements) result(this)
  implicit none
  
  type(ComplexSingleModeDisplacement), intent(in) :: displacements(:)
  type(ComplexModeDisplacement)                   :: this
  
  this%displacements = displacements
end function

! Return the number of modes along which the vector has displacements.
function size_ComplexModeDisplacement(input) result(output)
  implicit none
  
  type(ComplexModeDisplacement), intent(in) :: input
  integer                                   :: output
  
  output = size(input%displacements)
end function

! Returns a list of the q-points at which the displacement is non-zero.
function qpoints_ComplexModeDisplacement(this,complex_modes,qpoints) &
   & result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(ComplexMode),              intent(in) :: complex_modes(:)
  type(QpointData),               intent(in) :: qpoints(:)
  type(QpointData), allocatable              :: output(:)
  
  integer, allocatable :: qpoint_ids(:)
  
  integer :: i,j,ialloc
  
  allocate(qpoint_ids(size(this%displacements)), stat=ialloc); call err(ialloc)
  do i=1,size(this%displacements)
    j = first(complex_modes%id==this%displacements(i)%id)
    qpoint_ids(i) = complex_modes(j)%id
  enddo
  
  qpoint_ids = qpoint_ids(set(qpoint_ids))
  
  allocate(output(size(qpoint_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = qpoints(first(qpoints%id==qpoint_ids(i)))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_ComplexModeDisplacement(this,input)
  implicit none
  
  class(ComplexModeDisplacement), intent(out) :: this
  type(String),                   intent(in)  :: input(:)
  
  select type(this); type is(ComplexModeDisplacement)
    this = ComplexModeDisplacement(ComplexSingleModeDisplacement(input))
  end select
end subroutine

function write_ComplexModeDisplacement(this) result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(String), allocatable                  :: output(:)
  
  select type(this); type is(ComplexModeDisplacement)
    output = str(this%displacements)
  end select
end function

impure elemental function new_ComplexModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(ComplexModeDisplacement) :: this
  
  this = input
end function
end module
