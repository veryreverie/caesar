! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use real_mode_submodule
  use real_single_mode_displacement_submodule
  implicit none
  
  private
  
  public :: RealModeDisplacement
  public :: size
  
  type, extends(Stringsable) :: RealModeDisplacement
    type(RealSingleModeDisplacement), allocatable :: displacements(:)
  contains
    procedure, public :: modes   => modes_RealModeDisplacement
    procedure, public :: qpoints => qpoints_RealModeDisplacement
    
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface size
    module procedure size_RealModeDisplacement
  end interface
contains

! Return the number of modes along which the vector has displacements.
function size_RealModeDisplacement(input) result(output)
  implicit none
  
  type(RealModeDisplacement), intent(in) :: input
  integer                                :: output
  
  output = size(input%displacements)
end function

! Returns a list of the modes at which the displacement is non-zero.
function modes_RealModeDisplacement(this,real_modes) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(RealMode),              intent(in) :: real_modes(:)
  type(RealMode), allocatable             :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = real_modes(first(real_modes%id==this%displacements(i)%id))
  enddo
end function

! Returns a list of the q-points at which the displacement is non-zero.
function qpoints_RealModeDisplacement(this,real_modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(RealMode),              intent(in) :: real_modes(:)
  type(QPointData),            intent(in) :: qpoints(:)
  type(QPointData), allocatable           :: output(:)
  
  type(RealMode), allocatable :: modes(:)
  
  integer, allocatable :: qpoint_ids(:)
  
  integer :: i,j,ialloc
  
  ! List the q-point IDs of the modes in the displacement.
  modes = this%modes(real_modes)
  qpoint_ids = modes%qpoint_id
  
  ! De-duplicate the q-point IDs.
  qpoint_ids = qpoint_ids(set(qpoint_ids))
  
  ! List the q-points matching the IDs.
  allocate(output(size(qpoint_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i) = qpoints(first(qpoints%id==qpoint_ids(i)))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeDisplacement(this,input)
  implicit none
  
  class(RealModeDisplacement), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealModeDisplacement)
    allocate(this%displacements(size(input)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      this%displacements(i) = input(i)
    enddo
  end select
end subroutine

function write_RealModeDisplacement(this) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  integer :: i,ialloc
  
  select type(this); type is(RealModeDisplacement)
    allocate(output(size(this)), stat=ialloc); call err(ialloc)
    do i=1,size(this)
      output(i) = str(this%displacements(i))
    enddo
  end select
end function
end module
