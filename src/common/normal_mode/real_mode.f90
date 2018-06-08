! ======================================================================
! As complex_mode_module, but in real co-ordinates.
! ======================================================================
module real_mode_submodule
  use utils_module
  
  use structure_module
  
  use cartesian_vector_submodule
  implicit none
  
  private
  
  public :: RealMode
  
  ! A normal mode in real co-ordinates.
  type, extends(Stringsable) :: RealMode
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The mode in primitive cartesian co-ordinates.
    type(RealVector), allocatable :: primitive_vectors(:)
    
    ! The ID of the q-point corresponding to a complex mode related to this
    !    mode.
    ! N.B. real modes are not localised to a single q-point, but rather
    !   to a +/- pair of q-points. This id points to one of those points.
    integer :: qpoint_id
    
    ! The ID of the subspace of modes which are degenerate with this mode.
    integer :: subspace_id
  contains
    ! Return the mode in the cartesian co-ordinates of a given supercell.
    procedure, public :: cartesian_vector => cartesian_vector_RealMode
    
    ! I/O.
    procedure, public :: read  => read_RealMode
    procedure, public :: write => write_RealMode
  end type
  
  interface RealMode
    module procedure new_RealMode
    module procedure new_RealMode_StringArray
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_RealMode(id,paired_id,frequency,soft_mode,translational_mode, &
   & primitive_vectors,qpoint_id,subspace_id) result(this)
  implicit none
  
  integer,          intent(in) :: id
  integer,          intent(in) :: paired_id
  real(dp),         intent(in) :: frequency
  logical,          intent(in) :: soft_mode
  logical,          intent(in) :: translational_mode
  type(RealVector), intent(in) :: primitive_vectors(:)
  integer,          intent(in) :: qpoint_id
  integer,          intent(in) :: subspace_id
  type(RealMode)               :: this
  
  this%id                 = id
  this%paired_id          = paired_id
  this%frequency          = frequency
  this%soft_mode          = soft_mode
  this%translational_mode = translational_mode
  this%primitive_vectors  = primitive_vectors
  this%qpoint_id          = qpoint_id
  this%subspace_id        = subspace_id
end function

! ----------------------------------------------------------------------
! Return the mode in the cartesian co-ordinates of a given supercell.
! ----------------------------------------------------------------------
function cartesian_vector_RealMode(this,structure,qpoint) &
   & result(output)
  implicit none
  
  class(RealMode),     intent(in) :: this
  type(StructureData), intent(in) :: structure
  type(QpointData),    intent(in) :: qpoint
  type(CartesianVector)           :: output
  
  type(RealVector), allocatable :: vectors(:)
  
  type(AtomData)   :: atom
  type(RealVector) :: prim_vector
  type(IntVector)  :: atom_rvector
  real(dp)         :: qr
  
  integer :: i,ialloc
  
  if (qpoint%id/=this%qpoint_id) then
    call print_line(CODE_ERROR//': Mode and q-point incompatible.')
    call err()
  endif
  
  allocate(vectors(structure%no_atoms), stat=ialloc); call err(ialloc)
  do i=1,size(vectors)
    atom         = structure%atoms(i)
    prim_vector  = this%primitive_vectors(atom%prim_id())
    atom_rvector = structure%rvectors(atom%rvec_id())
    
    qr = dble(qpoint%qpoint*atom_rvector)
    
    if (this%id<=this%paired_id) then
      ! This is a cos(2 pi q.R) mode.
      vectors(i) = prim_vector * cos(2*PI*qr)
    else
      ! This is a sin(2 pi q.R) mode.
      vectors(i) = prim_vector * sin(2*PI*qr)
    endif
  enddo
  
  output = CartesianVector(vectors)
end function

! ----------------------------------------------------------------------
! I/O.
! ----------------------------------------------------------------------
subroutine read_RealMode(this,input)
  implicit none
  
  class(RealMode), intent(out) :: this
  type(String),    intent(in)  :: input(:)
  
  integer                       :: id
  integer                       :: paired_id
  real(dp)                      :: frequency
  logical                       :: soft_mode
  logical                       :: translational_mode
  type(RealVector), allocatable :: primitive_vectors(:)
  integer                       :: qpoint_id
  integer                       :: subspace_id
  
  type(String), allocatable :: line(:)
  
  integer :: no_atoms
  
  integer :: i,ialloc
  
  select type(this); type is(RealMode)
    ! Read the id of this mode.
    line = split_line(input(1))
    id = int(line(4))
    
    ! Read the id of this mode's pair.
    line = split_line(input(2))
    paired_id = int(line(6))
    
    ! Read frequency.
    line = split_line(input(3))
    frequency = dble(line(4))
    
    ! Read whether or not this mode is soft.
    line = split_line(input(4))
    soft_mode = lgcl(line(5))
    
    ! Read whether or not this mode is purely translational.
    line = split_line(input(5))
    translational_mode = lgcl(line(5))
    
    ! Read the q-qpoint id of this mode.
    line = split_line(input(6))
    qpoint_id = int(line(4))
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(7))
    subspace_id = int(line(4))
    
    ! Read in the vector associated with the mode.
    no_atoms = size(input)-8
    allocate( primitive_vectors(no_atoms), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms
      line = split_line(input(8+i))
      primitive_vectors(i) = dble(line)
    enddo
    
    this = RealMode( id,                 &
                   & paired_id,          &
                   & frequency,          &
                   & soft_mode,          &
                   & translational_mode, &
                   & primitive_vectors,  &
                   & qpoint_id,          &
                   & subspace_id)
  end select
end subroutine

function write_RealMode(this) result(output)
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String), allocatable   :: output(:)
  
  integer :: i,ialloc
  
  
  select type(this); type is(RealMode)
    allocate( output(8+size(this%primitive_vectors)), &
            & stat=ialloc); call err(ialloc)
    output(1) = 'Mode ID                   : '//this%id
    output(2) = 'ID of paired mode         : '//this%paired_id
    output(3) = 'Mode frequency            : '//this%frequency
    output(4) = 'Mode is soft              : '//this%soft_mode
    output(5) = 'Mode purely translational : '//this%translational_mode
    output(6) = 'q-point id                : '//this%qpoint_id
    output(7) = 'Degeneracy id             : '//this%subspace_id
    output(8) = 'Mode in cartesian primitive cell co-ordinates:'
    do i=1,size(this%primitive_vectors)
      output(8+i) = str(this%primitive_vectors(i))
    enddo
  end select
end function

impure elemental function new_RealMode_StringArray(input) &
   & result(this)
  implicit none
  
  type(StringArray), intent(in) :: input
  type(RealMode)                :: this
  
  this = input
end function
end module
