! ======================================================================
! Holds information about a degenerate subspace.
! ======================================================================
module degeneracy_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  use printable_module
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
    ! mode_ids(pairs_(i)) gives the ID of the paired mode to mode_ids(i).
    integer, allocatable, private :: pairs_(:)
  contains
    procedure, public :: modes => modes_DegenerateModes
    procedure, public :: qpoints => qpoints_DegenerateModes
    
    procedure, public :: str => str_DegenerateModes
  end type
  
  interface size
    module procedure size_DegenerateModes
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function process_degeneracies(modes,mode_qpoints) result(output)
  use normal_mode_module
  use qpoints_module
  use logic_module
  implicit none
  
  type(ComplexMode), intent(in)      :: modes(:)
  integer,           intent(in)      :: mode_qpoints(:)
  type(DegenerateModes), allocatable :: output(:)
  
  integer, allocatable :: degeneracy_ids(:)
  
  integer :: i,j,ialloc
  
  ! Make a list of degeneracy ids, with each id appearing exactly once.
  degeneracy_ids = modes%degeneracy_id
  degeneracy_ids = degeneracy_ids(set(degeneracy_ids))
  
  allocate(output(size(degeneracy_ids)), stat=ialloc); call err(ialloc)
  do i=1,size(output)
    output(i)%id = degeneracy_ids(i)
    output(i)%modes_ = filter(modes%degeneracy_id==output(i)%id)
    output(i)%mode_ids = modes(output(i)%modes_)%id
    output(i)%qpoints_ = mode_qpoints(output(i)%modes_)
    allocate(output(i)%pairs_(size(output(i))), stat=ialloc); call err(ialloc)
    do j=1,size(output(i))
      output(i)%pairs_(j) = &
         & first(output(i)%mode_ids==modes(output(i)%modes_(j))%paired_id)
    enddo
    
    ! Check pairs. Each mode should be paired with exactly one other.
    if (any(output(i)%pairs_==0)) then
      call print_line(CODE_ERROR//': Error pairing degenerate modes.')
      call err()
    endif
    do j=1,size(output(i))
      if (output(i)%pairs_(output(i)%pairs_(j))/=j) then
        call print_line(CODE_ERROR//': Error pairing degenerate modes.')
        call err()
      endif
    enddo
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
  use normal_mode_module
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
  use qpoints_module
  implicit none
  
  class(DegenerateModes), intent(in) :: this
  type(QpointData),       intent(in) :: qpoints(:)
  type(QpointData), allocatable      :: output(:)
  
  output = qpoints(this%qpoints_)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
recursive function str_DegenerateModes(this) result(output)
  implicit none
  
  class(DegenerateModes), intent(in) :: this
  type(String), allocatable          :: output(:)
  
  integer :: ialloc
  
  allocate(output(2), stat=ialloc); call err(ialloc)
  output(1) = 'Degeneracy ID : '//this%id
  output(2) = 'Mode IDs      : '//join(this%mode_ids)
end function
end module
