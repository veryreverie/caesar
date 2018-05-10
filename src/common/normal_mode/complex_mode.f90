! ======================================================================
! Harmonic normal modes, in complex co-ordinates.
! ======================================================================
module complex_mode_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: ComplexMode
  public :: operator(*)
  public :: l2_norm
  public :: conjg
  
  ! A normal mode in complex co-ordinates.
  type ComplexMode
    ! An id which is unique to each mode, and the unique id of the equivalent
    !    mode at the q-point q' s.t. q+q' is a G-vector.
    ! If 2q is a G-vector, then paired_id=id.
    integer :: id
    integer :: paired_id
    
    ! The frequency, and frequency-relevant information.
    real(dp) :: frequency
    logical  :: soft_mode          ! True if frequency < -epsilon.
    logical  :: translational_mode ! True if frequency=0 and at gamma.
    
    ! The displacements of atoms in the primitive cell.
    type(ComplexVector), allocatable :: primitive_displacements(:)
    
    ! An id which is shared between degenerate states, and different otherwise.
    integer :: degeneracy_id
  contains
    procedure, public :: write_file => write_file_ComplexMode
  end type
  
  interface ComplexMode
    module procedure read_file_ComplexMode
  end interface
  
  interface operator(*)
    module procedure dot_ComplexModes
  end interface
  
  interface l2_norm
    module procedure l2_norm_ComplexMode
  end interface
  
  interface conjg
    module procedure conjg_ComplexMode
  end interface
contains

! ----------------------------------------------------------------------
! ComplexMode procedures.
! ----------------------------------------------------------------------
subroutine write_file_ComplexMode(this,filename)
  implicit none
  
  class(ComplexMode), intent(in) :: this
  type(String),       intent(in) :: filename
  
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

function read_file_ComplexMode(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(ComplexMode)        :: this
  
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
    this%primitive_displacements(i) = cmplx(line)
  enddo
end function

! ----------------------------------------------------------------------
! Find the dot product between two complex modes.
! ----------------------------------------------------------------------
impure elemental function dot_ComplexModes(this,that) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  type(ComplexMode), intent(in) :: that
  complex(dp)                   :: output
  
  output = sum( this%primitive_displacements &
            & * that%primitive_displacements )
end function

impure elemental function l2_norm_ComplexMode(this) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: this
  real(dp)                      :: output
  
  output = sqrt(real(this*conjg(this)))
end function

! ----------------------------------------------------------------------
! Take the conjugate of the mode.
! This is equivalent to the transformation q -> -q.
! ----------------------------------------------------------------------
impure elemental function conjg_ComplexMode(input) result(output)
  implicit none
  
  type(ComplexMode), intent(in) :: input
  type(ComplexMode)             :: output
  
  output                         = input
  output%id                      = input%paired_id
  output%paired_id               = input%id
  output%primitive_displacements = conjg(input%primitive_displacements)
end function
end module
