! ======================================================================
! Harmonic normal modes.
! ======================================================================
module normal_mode_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  type NormalMode
    real(dp)                         :: frequency
    type(ComplexVector), allocatable :: displacements(:)
    logical                          :: soft_mode
    logical                          :: translational_mode
  end type
  
  type ModeVector
    real(dp), allocatable :: vector(:)
  contains
    generic,   public  :: operator(+) => add_ModeVector_ModeVector
    procedure, private ::                add_ModeVector_ModeVector
    generic,   public  :: operator(-) => subtract_ModeVector_ModeVector
    procedure, private ::                subtract_ModeVector_ModeVector
  end type
contains

function add_ModeVector_ModeVector(a,b) result(output)
  implicit none
  
  class(ModeVector), intent(in) :: a
  class(ModeVector), intent(in) :: b
  type(ModeVector)              :: output
  
  output%vector = a%vector + b%vector
end function

function subtract_ModeVector_ModeVector(a,b) result(output)
  implicit none
  
  class(ModeVector), intent(in) :: a
  class(ModeVector), intent(in) :: b
  type(ModeVector)              :: output
  
  output%vector = a%vector - b%vector
end function

function read_normal_mode_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(NormalMode)         :: this
  
  type(String), allocatable :: mode_file(:)
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  
  integer :: i,ialloc
  
  ! Read in mode file.
  mode_file = read_lines(filename)
  
  ! Allocate space.
  no_atoms = size(mode_file)-4
  allocate(this%displacements(no_atoms), stat=ialloc); call err(ialloc)
  
  ! Read frequency.
  line = split(mode_file(1))
  this%frequency = dble(line(2))
  
  ! Read whether or not this mode is soft.
  line = split(mode_file(2))
  this%soft_mode = lgcl(line(3))
  
  ! Read whether or not this mode is purely translational.
  line = split(mode_file(3))
  this%translational_mode = lgcl(line(3))
  
  ! Read in the displacement associated with the mode.
  do i=1,no_atoms
    line = split(mode_file(4+i))
    this%displacements(i) = cmplx(line)
  enddo
end function

subroutine write_normal_mode_file(this,filename)
  implicit none
  
  type(NormalMode), intent(in) :: this
  type(String),     intent(in) :: filename
  
  integer :: mode_file
  integer :: i
  
  mode_file = open_write_file(filename)
  call print_line(mode_file, 'Frequency:          '//this%frequency)
  call print_line(mode_file, 'Soft mode:          '//this%soft_mode)
  call print_line(mode_file, 'Translational mode: '//this%translational_mode)
  call print_line(mode_file, 'Displacements:')
  do i=1,size(this%displacements)
    call print_line(mode_file, this%displacements(i))
  enddo
  close(mode_file)
end subroutine

! ----------------------------------------------------------------------
! Converts a vector in normal mode co-ordinates to cartesian co-ordinates.
! ----------------------------------------------------------------------
function normal_mode_to_cartesian(input,modes,qpoint,supercell) result(output)
  use constants_module, only : pi
  use qpoints_module
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(ModeVector),    intent(in) :: input
  type(NormalMode),    intent(in) :: modes(:)
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: supercell
  type(RealVector), allocatable   :: output(:)
  
  type(RealVector)       :: direction
  
  ! Working variables.
  real(dp)    :: qr
  complex(dp) :: exp_iqr
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  allocate(output(size(supercell%atoms)), stat=ialloc); call err(ialloc)
  do i=1,size(supercell%atoms)
    output(i) = dble(int(zeroes(3)))
    
    ! Calculate 2*pi*q.R and exp(2*pi*i*q.R).
    qr = 2 * pi * qpoint%qpoint * supercell%rvectors(supercell%atom_to_rvec(i))
    exp_iqr = cmplx(cos(qr), sin(qr), dp)
    
    ! Calculate displacements in cartesian co-ordinates.
    do j=1,size(modes)
      ! Calculate the direction of the displacement.
      direction = real( modes(j)%displacements(supercell%atom_to_prim(i)) &
                    & * exp_iqr)
      ! Calculate the displacement.
      output(i) = output(i) + direction * input%vector(j)
    enddo
  enddo
end function

! ----------------------------------------------------------------------
! Converts a vector in cartesion co-ordinates to normal mode co-ordinates.
! ----------------------------------------------------------------------
function cartesian_to_normal_mode(input,modes,qpoint,supercell) result(output)
  use constants_module, only : pi
  use qpoints_module
  use structure_module
  use linear_algebra_module
  implicit none
  
  type(RealVector),    intent(in) :: input(:)
  type(NormalMode),    intent(in) :: modes(:)
  type(QpointData),    intent(in) :: qpoint
  type(StructureData), intent(in) :: supercell
  type(ModeVector)                :: output
  
  type(RealVector)       :: direction
  
  ! Working variables.
  real(dp)    :: qr
  complex(dp) :: exp_iqr
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  allocate(output%vector(size(modes)), stat=ialloc); call err(ialloc)
  output%vector = 0
  do i=1,size(supercell%atoms)
    ! Calculate 2*pi*q.R and exp(2*pi*i*q.R).
    qr = 2 * pi * qpoint%qpoint * supercell%rvectors(supercell%atom_to_rvec(i))
    exp_iqr = cmplx(cos(qr), sin(qr), dp)
    
    ! Calculate displacements in normal-mode co-ordinates.
    do j=1,size(modes)
      ! Calculate the direction of the displacement.
      direction = real( modes(j)%displacements(supercell%atom_to_prim(i)) &
                    & * exp_iqr)
      
      ! Calculate the projection of the vector onto the normal mode.
      output%vector(j) = output%vector(j) &
                     & + direction * input(i)
    enddo
  enddo
end function
end module
