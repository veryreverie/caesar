! ======================================================================
! Holds information about a degenerate subspace.
! ======================================================================
module caesar_degenerate_subspace_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  use caesar_real_mode_module
  implicit none
  
  private
  
  public :: DegenerateSubspace
  public :: process_degeneracies
  public :: size
  
  type, extends(Stringsable) :: DegenerateSubspace
    ! The id of the degeneracy.
    integer, public :: id
    
    ! The harmonic frequency of the modes in the degenerate subspace.
    real(dp), public :: frequency
    
    ! The IDs of the modes in the degenerate subspace.
    integer, allocatable, public :: mode_ids(:)
    
    ! paired_ids(i) is the ID of the mode paired to mode_ids(i).
    integer, allocatable, public :: paired_ids(:)
  contains
    generic,   public  :: modes =>                               &
                        & modes_DegenerateSubspace_ComplexModes, &
                        & modes_DegenerateSubspace_RealModes
    procedure, private :: modes_DegenerateSubspace_ComplexModes
    procedure, private :: modes_DegenerateSubspace_RealModes
    
    generic,   public  :: qpoints =>                               &
                        & qpoints_DegenerateSubspace_ComplexModes, &
                        & qpoints_DegenerateSubspace_RealModes
    procedure, private :: qpoints_DegenerateSubspace_ComplexModes
    procedure, private :: qpoints_DegenerateSubspace_RealModes
    
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
function new_DegenerateSubspace(id,frequency,mode_ids,paired_ids) result(this)
  implicit none
  
  integer,  intent(in)     :: id
  real(dp), intent(in)     :: frequency
  integer,  intent(in)     :: mode_ids(:)
  integer,  intent(in)     :: paired_ids(:)
  type(DegenerateSubspace) :: this
  
  if (size(mode_ids)/=size(paired_ids)) then
    call print_line(CODE_ERROR//': mode IDs and paired IDs do not match.')
    call err()
  endif
  
  this%id         = id
  this%frequency  = frequency
  this%mode_ids   = mode_ids
  this%paired_ids = paired_ids
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
    
    output(i) = DegenerateSubspace( id         = subspace_ids(i),             &
                                  & frequency  = subspace_modes(1)%frequency, &
                                  & mode_ids   = subspace_modes%id,           &
                                  & paired_ids = subspace_modes%paired_id     )
  enddo
end function

! ----------------------------------------------------------------------
! Returns the degenerate modes.
! ----------------------------------------------------------------------
function modes_DegenerateSubspace_ComplexModes(this,modes) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(ComplexMode),         intent(in) :: modes(:)
  type(ComplexMode), allocatable        :: output(:)
  
  integer :: i,j,ialloc
  
  !output = [( modes(first(modes%id==this%mode_ids(i))), i=1, size(this) )]
  ! WORKAROUND: To avoid a memory leak in ifort 19.1.0.166
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    j = first(modes%id==this%mode_ids(i))
    output(i) = modes(j)
  enddo
end function

function modes_DegenerateSubspace_RealModes(this,modes) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(RealMode),            intent(in) :: modes(:)
  type(RealMode), allocatable           :: output(:)
  
  integer :: i,j,ialloc
  
  !output = [( modes(first(modes%id==this%mode_ids(i))), i=1, size(this) )]
  ! WORKAROUND: To avoid a memory leak in ifort 19.1.0.166
  allocate(output(size(this)), stat=ialloc); call err(ialloc)
  do i=1,size(this)
    j = first(modes%id==this%mode_ids(i))
    output(i) = modes(j)
  enddo
end function

! ----------------------------------------------------------------------
! Returns the q-points corresponding to the degenerate modes.
! ----------------------------------------------------------------------
function qpoints_DegenerateSubspace_ComplexModes(this,modes,qpoints) &
   & result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(ComplexMode),         intent(in) :: modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(QpointData), allocatable         :: output(:)
  
  type(ComplexMode), allocatable :: subspace_modes(:)
  
  integer :: i
  
  subspace_modes = this%modes(modes)
  output = [( qpoints(first(qpoints%id==subspace_modes(i)%qpoint_id)), &
            & i=1,                                                     &
            & size(subspace_modes)                                     )]
end function

function qpoints_DegenerateSubspace_RealModes(this,modes,qpoints) &
   & result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(RealMode),            intent(in) :: modes(:)
  type(QpointData),          intent(in) :: qpoints(:)
  type(QpointData), allocatable         :: output(:)
  
  type(RealMode), allocatable :: subspace_modes(:)
  
  integer :: i
  
  subspace_modes = this%modes(modes)
  output = [( qpoints(first(qpoints%id==subspace_modes(i)%qpoint_id_plus)), &
            & i=1,                                                          &
            & size(subspace_modes)                                          )]
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_DegenerateSubspace(this,input)
  implicit none
  
  class(DegenerateSubspace), intent(out) :: this
  type(String),              intent(in)  :: input(:)
  
  integer              :: id
  real(dp)             :: frequency
  integer, allocatable :: mode_ids(:)
  integer, allocatable :: paired_ids(:)
  
  type(String), allocatable :: line(:)
  
  select type(this); type is(DegenerateSubspace)
    line = split_line(input(1))
    id = int(line(5))
    
    line = split_line(input(2))
    frequency = dble(line(5))
    
    line = split_line(input(3))
    mode_ids = int(line(5:))
    
    line = split_line(input(4))
    paired_ids = int(line(5:))
    
    this = DegenerateSubspace(id, frequency, mode_ids, paired_ids)
  class default
    call err()
  end select
end subroutine

function write_DegenerateSubspace(this) result(output)
  implicit none
  
  class(DegenerateSubspace), intent(in) :: this
  type(String), allocatable             :: output(:)
  
  select type(this); type is(DegenerateSubspace)
    output = [ 'Degenerate subspace ID : '//this%id,        &
             & 'Frequency of modes     : '//this%frequency, &
             & 'Degenerate mode IDs    : '//this%mode_ids,  &
             & 'Paired mode IDS        : '//this%paired_ids ]
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
