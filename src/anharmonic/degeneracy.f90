! ======================================================================
! Holds information about a degenerate subspace.
! ======================================================================
module degeneracy_module
  use common_module
  implicit none
  
  private
  
  public :: DegenerateModes
  public :: process_degeneracies
  public :: size
  
  type, extends(printable) :: DegenerateModes
    ! --------------------------------------------------
    ! Public variables, corresponding to class ids.
    ! --------------------------------------------------
    ! The id of the degeneracy.
    integer, public :: id
    
    ! The ids of the modes in the degenerate subspace.
    integer, allocatable, public :: mode_ids(:)
    
    ! --------------------------------------------------
    ! Private variables, corresponding to the positions of
    !    objects within lists.
    ! --------------------------------------------------
    ! normal_modes(modes_(i)) is the i'th degenerate mode.
    ! This mode has ID mode_ids(i).
    integer, allocatable, private :: modes_(:)
    ! qpoints(qpoints_(i)) is the q-point of the i'th degenerate mode.
    integer, allocatable, private :: qpoints_(:)
  contains
    procedure, public :: modes => modes_DegenerateModes
    procedure, public :: qpoints => qpoints_DegenerateModes
    
    procedure, public :: to_String => to_String_DegenerateModes
  end type
  
  interface size
    module procedure size_DegenerateModes
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function process_degeneracies(modes,mode_qpoints) result(output)
  implicit none
  
  type(ComplexMode), intent(in)      :: modes(:)
  integer,           intent(in)      :: mode_qpoints(:)
  type(DegenerateModes), allocatable :: output(:)
  
  integer, allocatable :: degeneracy_ids(:)
  
  integer :: i,j,ialloc
  
  ! Make a list of degeneracy ids.
  degeneracy_ids = modes%degeneracy_id
  ! Remove the purely translational modes.
  degeneracy_ids = degeneracy_ids(filter(.not.modes%translational_mode))
  ! De-duplicate the list, so that each id appears exactly once.
  degeneracy_ids = degeneracy_ids(set(degeneracy_ids))
  
  ! Generate subspaces.
  allocate(output(size(degeneracy_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%id = degeneracy_ids(i)
    output(i)%modes_ = filter(modes%degeneracy_id==output(i)%id)
    output(i)%mode_ids = modes(output(i)%modes_)%id
    output(i)%qpoints_ = mode_qpoints(output(i)%modes_)
  enddo
end function

! ----------------------------------------------------------------------
! Returns the number of degenerate modes.
! ----------------------------------------------------------------------
function size_DegenerateModes(input) result(output)
  implicit none
  
  type(DegenerateModes), intent(in) :: input
  integer                           :: output
  
  output = size(input%mode_ids)
end function

! ----------------------------------------------------------------------
! Returns the degenerate modes.
! ----------------------------------------------------------------------
function modes_DegenerateModes(this,modes) result(output)
  implicit none
  
  class(DegenerateModes), intent(in) :: this
  type(ComplexMode),      intent(in) :: modes(:)
  type(ComplexMode), allocatable     :: output(:)
  
  output = modes(this%modes_)
end function

! ----------------------------------------------------------------------
! Returns the q-points corresponding to the degenerate modes.
! ----------------------------------------------------------------------
function qpoints_DegenerateModes(this,qpoints) result(output)
  implicit none
  
  class(DegenerateModes), intent(in) :: this
  type(QpointData),       intent(in) :: qpoints(:)
  type(QpointData), allocatable      :: output(:)
  
  output = qpoints(this%qpoints_)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
recursive function to_String_DegenerateModes(this) result(output)
  implicit none
  
  class(DegenerateModes), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  integer :: ialloc
  
  allocate(output(2), stat=ialloc); call err(ialloc)
  output(1) = 'Degeneracy ID : '//this%id
  output(2) = 'Mode IDs      : '//join(this%mode_ids)
end function
end module
