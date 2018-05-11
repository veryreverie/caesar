! ======================================================================
! As complex_mode_module, but in real co-ordinates.
! ======================================================================
module real_mode_submodule
  use utils_module
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
    
    ! The displacements of atoms in the primitive cell.
    type(RealVector), allocatable :: primitive_displacements(:)
    
    ! An id which is shared between degenerate states, and different otherwise.
    integer :: degeneracy_id
  contains
    procedure, public :: read  => read_RealMode
    procedure, public :: write => write_RealMode
  end type
  
  interface RealMode
    module procedure new_RealMode
  end interface
contains

! ----------------------------------------------------------------------
! Constructor.
! ----------------------------------------------------------------------
function new_RealMode(id,paired_id,frequency,soft_mode,translational_mode, &
   & primitive_displacements, degeneracy_id) result(this)
  implicit none
  
  integer,          intent(in) :: id
  integer,          intent(in) :: paired_id
  real(dp),         intent(in) :: frequency
  logical,          intent(in) :: soft_mode
  logical,          intent(in) :: translational_mode
  type(RealVector), intent(in) :: primitive_displacements(:)
  integer,          intent(in) :: degeneracy_id
  type(RealMode)               :: this
  
  this%id = id
  this%paired_id = paired_id
  this%frequency = frequency
  this%soft_mode = soft_mode
  this%translational_mode = translational_mode
  this%primitive_displacements = primitive_displacements
  this%degeneracy_id = degeneracy_id
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
  type(RealVector), allocatable :: primitive_displacements(:)
  integer                       :: degeneracy_id
  
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
    
    ! Read the degeneracy id of this mode.
    line = split_line(input(6))
    degeneracy_id = int(line(4))
    
    ! Read in the displacement associated with the mode.
    no_atoms = size(input)-7
    allocate( primitive_displacements(no_atoms), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_atoms
      line = split_line(input(7+i))
      primitive_displacements(i) = dble(line)
    enddo
    
    this = RealMode( id,                      &
                      & paired_id,               &
                      & frequency,               &
                      & soft_mode,               &
                      & translational_mode,      &
                      & primitive_displacements, &
                      & degeneracy_id)
  end select
end subroutine

function write_RealMode(this) result(output)
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String), allocatable   :: output(:)
  
  integer :: i,ialloc
  
  
  select type(this); type is(RealMode)
    allocate( output(7+size(this%primitive_displacements)), &
            & stat=ialloc); call err(ialloc)
    output(1) = 'Mode ID                   : '//this%id
    output(2) = 'ID of paired mode         : '//this%paired_id
    output(3) = 'Mode frequency            : '//this%frequency
    output(4) = 'Mode is soft              : '//this%soft_mode
    output(5) = 'Mode purely translational : '//this%translational_mode
    output(6) = 'Degeneracy id             : '//this%degeneracy_id
    output(7) = 'Displacements in primitive cell:'
    do i=1,size(this%primitive_displacements)
      output(7+i) = str(this%primitive_displacements(i))
    enddo
  end select
end function
end module
