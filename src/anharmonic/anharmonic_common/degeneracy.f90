! ======================================================================
! Holds information about a degenerate subspace.
! ======================================================================
module degeneracy_module
  use common_module
  implicit none
  
  private
  
  public :: DegenerateSubspace
  public :: process_degeneracies
  public :: size
  
  type, extends(Stringsable) :: DegenerateSubspace
    ! The id of the degeneracy.
    integer, public :: id
    
    ! The ids of the modes in the degenerate subspace.
    integer, allocatable, public :: mode_ids(:)
  contains
    procedure, public :: modes => modes_DegenerateSubspace
    procedure, public :: qpoints => qpoints_DegenerateSubspace
    
    ! I/O.
    procedure, public :: read  => read_DegenerateSubspace
    procedure, public :: write => write_DegenerateSubspace
  end type
  
  interface DegenerateSubspace
    module procedure new_DegenerateSubspace
    module procedure new_DegenerateSubspace_Strings
    module procedure new_DegenerateSubspace_StringArray
  end interface
  
  interface size
    module procedure size_DegenerateSubspace
  end interface
contains

! ----------------------------------------------------------------------
! Basic functionality: constructor and size() function.
! ----------------------------------------------------------------------
function new_DegenerateSubspace(id,mode_ids) result(this)
  implicit none
  
  integer, intent(in)      :: id
  integer, intent(in)      :: mode_ids(:)
  type(DegenerateSubspace) :: this
  
  this%id       = id
  this%mode_ids = mode_ids
end function

function size_DegenerateSubspace(input) result(output)
  implicit none
  
  type(DegenerateSubspace), intent(in) :: input
  integer                              :: output
  
  output = size(input%mode_ids)
end function

! ----------------------------------------------------------------------
! Construct the degenerate subspaces.
! ----------------------------------------------------------------------
function process_degeneracies(modes) result(output)
  implicit none
  
  type(ComplexMode), intent(in)         :: modes(:)
  type(DegenerateSubspace), allocatable :: output(:)
  
  integer,           allocatable :: subspace_ids(:)
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer :: i,ialloc
  
  ! Make a list of degeneracy ids.
  subspace_ids = modes%subspace_id
  ! Remove the purely translational modes.
  subspace_ids = subspace_ids(filter(.not.modes%translational_mode))
  ! De-duplicate the list, so that each id appears exactly once.
  subspace_ids = subspace_ids(set(subspace_ids))
  
  ! Generate subspaces.
  allocate(output(size(subspace_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    subspace_modes = modes(filter(modes%subspace_id==subspace_ids(i)))
    output(i)%id = subspace_ids(i)
    output(i)%mode_ids = subspace_modes%id
  enddo
end function

! ----------------------------------------------------------------------
! Returns the degenerate modes.
! ----------------------------------------------------------------------
function modes_DegenerateSubspace(this,modes) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(ComplexMode),         intent(in) :: modes(:)
  type(ComplexMode), allocatable        :: output(:)
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    output(i) = modes(first(modes%id==this%mode_ids(i)))
  enddo
end function

! ----------------------------------------------------------------------
! Returns the q-points corresponding to the degenerate modes.
! ----------------------------------------------------------------------
function qpoints_DegenerateSubspace(this,modes,qpoints) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(ComplexMode),         intent(in) :: modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(QpointData), allocatable         :: output(:)
  
  type(ComplexMode) :: mode
  
  integer :: i,ialloc
  
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    mode = modes(first(modes%id==this%mode_ids(i)))
    output(i) = qpoints(first(qpoints%id==mode%qpoint_id))
  enddo
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_DegenerateSubspace(this,input)
  implicit none
  
  class(DegenerateSubspace), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  integer              :: id
  integer, allocatable :: mode_ids(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(DegenerateSubspace)
    line = split_line(input(1))
    id = int(line(5))
    
    line = split_line(input(2))
    mode_ids = int(line(5:))
    
    this = DegenerateSubspace(id,mode_ids)
  class default
    call err()
  end select
end subroutine

function write_DegenerateSubspace(this) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(DegenerateSubspace)
    output = [ 'Degenerate subspace ID : '//this%id, &
             & 'Degenerate mode IDs    : '//this%mode_ids ]
  class default
    call err()
  end select
end function

function new_DegenerateSubspace_Strings(input) result(this)
  implicit none
  
  type(String), intent(in) :: input(:)
  type(DegenerateSubspace) :: this
  
  call this%read(input)
end function

impure elemental function new_DegenerateSubspace_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(DegenerateSubspace)      :: this
  
  this = DegenerateSubspace(str(input))
end function
end module
