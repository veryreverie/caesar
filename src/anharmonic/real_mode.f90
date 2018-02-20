! ======================================================================
! As normal_mode_module, but in real co-ordinates.
! ======================================================================
module real_mode_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  private
  
  public :: RealMode
  public :: complex_to_real
  
  ! A normal mode in real co-ordinates.
  type RealMode
    ! Whether or not 2q=G. If true, there is only the cosine mode.
    !    If not, there is a cosine mode and a sine mode.
    logical :: at_paired_qpoint
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The displacements of atoms in the primitive cell.
    type(RealVector), allocatable :: primitive_displacements(:)
  contains
    procedure, public :: write_file => write_file_RealMode
  end type
  
  interface RealMode
    module procedure read_file_RealMode
  end interface
  
  ! Conversions between complex and real co-ordinates.
  interface complex_to_real
    module procedure complex_to_real_Mode
  end interface
contains
! ----------------------------------------------------------------------
! RealMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_RealMode(this,filename)
  use ofile_module
  implicit none
  
  class(RealMode), intent(in) :: this
  type(String),    intent(in) :: filename
  
  type(OFile) :: mode_file
  
  integer :: i
  
  mode_file = filename
  call mode_file%print_line('Mode at paired q-point: '//this%at_paired_qpoint)
  call mode_file%print_line('Frequency:              '//this%frequency)
  call mode_file%print_line('Soft mode:              '//this%soft_mode)
  call mode_file%print_line('Translational mode:     '//this%translational_mode)
  call mode_file%print_line('Displacements in primitive cell:')
  do i=1,size(this%primitive_displacements)
    call mode_file%print_line(this%primitive_displacements(i))
  enddo
end subroutine

function read_file_RealMode(filename) result(this)
  use ifile_module
  implicit none
  
  type(String), intent(in) :: filename
  type(RealMode)           :: this
  
  type(IFile)               :: mode_file
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = filename
  
  ! Read kind.
  line = split(mode_file%line(1))
  this%at_paired_qpoint = lgcl(line(5))
  
  ! Read frequency.
  line = split(mode_file%line(2))
  this%frequency = dble(line(2))
  
  ! Read whether or not this mode is soft.
  line = split(mode_file%line(3))
  this%soft_mode = lgcl(line(3))
  
  ! Read whether or not this mode is purely translational.
  line = split(mode_file%line(4))
  this%translational_mode = lgcl(line(3))
  
  ! Read in the displacement associated with the mode.
  no_atoms = size(mode_file)-5
  allocate( this%primitive_displacements(no_atoms), &
          & stat=ialloc); call err(ialloc)
  do i=1,no_atoms
    line = split(mode_file%line(5+i))
    this%primitive_displacements(i) = dble(line)
  enddo
end function

! Returns a cos mode and a sin mode if 2q/=G.
function complex_to_real_Mode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(RealMode), allocatable   :: output(:)
  
  integer :: no_atoms_prim
  
  integer :: i,ialloc
  
  ! Copy over common information.
  if (input%at_paired_qpoint) then
    allocate(output(1), stat=ialloc); call err(ialloc)
  else
    allocate(output(2), stat=ialloc); call err(ialloc)
  endif
  
  output%at_paired_qpoint   = input%at_paired_qpoint
  output%frequency          = input%frequency
  output%soft_mode          = input%soft_mode
  output%translational_mode = input%translational_mode
  
  no_atoms_prim = size(input%primitive_displacements)
  
  ! Convert displacements.
  if (input%at_paired_qpoint) then
    ! x = x.
    output(1)%primitive_displacements = real(input%primitive_displacements)
  else
    ! c = (x+ + x-)/sqrt(2) = real(x+) * sqrt(2).
    output(1)%primitive_displacements = real(input%primitive_displacements) &
                                    & * sqrt(2.0_dp)
    ! s = (x+ - x-)/(sqrt(2)i) = imag(x+) * sqrt(2).
    output(2)%primitive_displacements = aimag(input%primitive_displacements) &
                                    & * sqrt(2.0_dp)
  endif
end function
end module
