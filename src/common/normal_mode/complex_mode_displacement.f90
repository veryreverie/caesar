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
  
  type, extends(Stringsable) :: ComplexModeDisplacement
    type(ComplexSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: qpoints => qpoints_ComplexModeDisplacement
    
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

! Returns a list of the q-points at which the displacement is non-zero.
function qpoints_ComplexModeDisplacement(this,complex_modes,qpoints) &
   & result(output)
  implicit none
  
  class(ComplexModeDisplacement), intent(in) :: this
  type(ComplexMode),              intent(in) :: complex_modes(:)
  type(QPointData),               intent(in) :: qpoints(:)
  type(QPointData), allocatable              :: output(:)
  
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
