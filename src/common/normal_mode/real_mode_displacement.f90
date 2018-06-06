! ======================================================================
! A displacement in real mode co-ordinates.
! ======================================================================
module real_mode_displacement_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_displacement_submodule
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
    
    procedure, public :: cartesian_displacement => &
       & cartesian_displacement_RealModeDisplacement
    
    procedure, public :: read  => read_RealModeDisplacement
    procedure, public :: write => write_RealModeDisplacement
  end type
  
  interface RealModeDisplacement
    module procedure new_RealModeDisplacement
    module procedure new_RealModeDisplacement_StringArray
  end interface
  
  interface size
    module procedure size_RealModeDisplacement
  end interface
contains

! Constructor.
function new_RealModeDisplacement(displacements) result(this)
  implicit none
  
  type(RealSingleModeDisplacement), intent(in) :: displacements(:)
  type(RealModeDisplacement)                   :: this
  
  this%displacements = displacements
end function

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
  type(QpointData),            intent(in) :: qpoints(:)
  type(QpointData), allocatable           :: output(:)
  
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

! Returns the displacement in cartesian co-ordinates.
function cartesian_displacement_RealModeDisplacement(this,structure, &
   & real_modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: real_modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(CartesianDisplacement)             :: output
  
  type(RealMode),   allocatable :: modes(:)
  type(QpointData), allocatable :: mode_qpoints(:)
  type(IntVector),  allocatable :: mode_rvectors(:)
  
  type(CartesianDisplacement), allocatable :: mode_displacements(:)
  
  integer :: i,j,ialloc
  
  ! List the modes which have non-zero displacement,
  !   and the q-points and R-vectors associated with those modes.
  allocate( modes(size(this)),         &
          & mode_qpoints(size(this)),  &
          & mode_rvectors(size(this)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    j = first(real_modes%id==this%displacements(i)%id)
    modes(i) = real_modes(j)
    mode_qpoints(i) = qpoints(first(qpoints%id==modes(i)%qpoint_id))
  enddo
  
  ! Calculate the cartesian displacement due to each mode.
  allocate(mode_displacements(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    mode_displacements(i) =                                       &
       & this%displacements(i)%cartesian_displacement( modes(i),  &
       &                                               structure, &
       &                                               mode_qpoints(i))
  enddo
  
  ! Add together the contributions from all modes.
  output = sum(mode_displacements)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealModeDisplacement(this,input)
  implicit none
  
  class(RealModeDisplacement), intent(out) :: this
  type(String),                intent(in)  :: input(:)
  
  select type(this); type is(RealModeDisplacement)
    this = RealModeDisplacement(RealSingleModeDisplacement(input))
  end select
end subroutine

function write_RealModeDisplacement(this) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(String), allocatable               :: output(:)
  
  select type(this); type is(RealModeDisplacement)
    output = str(this%displacements)
  end select
end function

impure elemental function new_RealModeDisplacement_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealModeDisplacement)    :: this
  
  this = input
end function
end module
