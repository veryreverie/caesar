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
    module procedure new_RealModeDisplacement_CartesianDisplacement
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

! ----------------------------------------------------------------------
! Conversions between CartesianDisplacement and RealModeDisplacement.
! ----------------------------------------------------------------------
! Returns the displacement in cartesian co-ordinates.
function cartesian_displacement_RealModeDisplacement(this,structure, &
   & modes,qpoints) result(output)
  implicit none
  
  class(RealModeDisplacement), intent(in) :: this
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(CartesianDisplacement)             :: output
  
  type(RealMode)   :: mode
  type(QpointData) :: qpoint
  
  type(CartesianDisplacement), allocatable :: displacements(:)
  
  integer :: i,ialloc
  
  ! Calculate the cartesian displacement due to each mode.
  allocate(displacements(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    ! Find the mode and q-point associated with displacement i.
    mode = modes(first(modes%id==this%displacements(i)%id))
    qpoint = qpoints(first(qpoints%id==mode%qpoint_id))
    
    ! Calculate the displacement from mode i.
    displacements(i) =                                            &
       & this%displacements(i)%cartesian_displacement( mode,      &
       &                                               structure, &
       &                                               qpoint)
  enddo
  
  ! Add together the contributions from all modes.
  output = sum(displacements)
end function

! Converts a CartesianDisplacement to a RealModeDisplacement.
function new_RealModeDisplacement_CartesianDisplacement(displacement, &
   & structure,modes,qpoints) result(this)
  implicit none
  
  type(CartesianDisplacement), intent(in) :: displacement
  type(StructureData),         intent(in) :: structure
  type(RealMode),              intent(in) :: modes(:)
  type(QpointData),            intent(in) :: qpoints(:)
  type(RealModeDisplacement)              :: this
  
  type(RealSingleModeDisplacement), allocatable :: displacements(:)
  type(QpointData)                              :: qpoint
  
  integer :: i,ialloc
  
  allocate(displacements(size(modes)), stat=ialloc); call err(ialloc)
  do i=1,size(modes)
    qpoint = qpoints(first(qpoints%id==modes(i)%qpoint_id))
    displacements(i) = RealSingleModeDisplacement( modes(i),     &
                                                 & displacement, &
                                                 & structure,    &
                                                 & qpoint)
  enddo
  
  this = RealModeDisplacement(displacements)
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
