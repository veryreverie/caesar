! ======================================================================
! As normal_mode_module, but in real co-ordinates.
! ======================================================================
module real_mode_submodule
  use utils_module
  implicit none
  
  private
  
  public :: RealMode
  
  ! A normal mode in real co-ordinates.
  type RealMode
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
    procedure, public :: write_file => write_file_RealMode
  end type
  
  interface RealMode
    module procedure read_file_RealMode
  end interface
contains

! ----------------------------------------------------------------------
! RealMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_RealMode(this,filename)
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String),    intent(in) :: filename
  
  type(OFile) :: mode_file
  
  integer :: i
  
  mode_file = OFile(filename)
  call mode_file%print_line('Mode ID                   : '//this%id)
  call mode_file%print_line('ID of paired mode         : '//this%paired_id)
  call mode_file%print_line('Mode frequency            : '//this%frequency)
  call mode_file%print_line('Mode is soft              : '//this%soft_mode)
  call mode_file%print_line('Mode purely translational : '// &
     & this%translational_mode)
  call mode_file%print_line('Degeneracy id             : '//this%degeneracy_id)
  call mode_file%print_line('Displacements in primitive cell:')
  do i=1,size(this%primitive_displacements)
    call mode_file%print_line(this%primitive_displacements(i))
  enddo
end subroutine

function read_file_RealMode(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(RealMode)           :: this
  
  type(IFile)               :: mode_file
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = IFile(filename)
  
  ! Read the id of this mode.
  line = split_line(mode_file%line(1))
  this%id = int(line(4))
  
  ! Read the id of this mode's pair.
  line = split_line(mode_file%line(2))
  this%paired_id = int(line(6))
  
  ! Read frequency.
  line = split_line(mode_file%line(3))
  this%frequency = dble(line(4))
  
  ! Read whether or not this mode is soft.
  line = split_line(mode_file%line(4))
  this%soft_mode = lgcl(line(5))
  
  ! Read whether or not this mode is purely translational.
  line = split_line(mode_file%line(5))
  this%translational_mode = lgcl(line(5))
  
  ! Read the degeneracy id of this mode.
  line = split_line(mode_file%line(6))
  this%degeneracy_id = int(line(4))
  
  ! Read in the displacement associated with the mode.
  no_atoms = size(mode_file)-7
  allocate( this%primitive_displacements(no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    line = split_line(mode_file%line(7+i))
    this%primitive_displacements(i) = dble(line)
  enddo
end function
end module
