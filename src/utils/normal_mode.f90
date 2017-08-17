! ======================================================================
! Harmonic normal modes.
! ======================================================================
module normal_mode_module
  use constants_module, only : dp
  use string_module
  use io_module
  use linear_algebra_module
  
  type NormalMode
    logical                          :: soft_mode
    real(dp)                         :: frequency
    type(ComplexVector), allocatable :: displacements(:)
  end type
  
  type ModeVector
    real(dp), allocatable :: vector(:)
  end type
  
  interface new
    module procedure new_NormalMode
  end interface
contains

subroutine new_NormalMode(this,no_atoms)
  implicit none
  
  type(NormalMode), intent(out) :: this
  integer         , intent(in)  :: no_atoms
  
  integer :: ialloc
  
  allocate(this%displacements(no_atoms), stat=ialloc); call err(ialloc)
end subroutine

function read_normal_mode_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(NormalMode)         :: this
  
  type(String), allocatable :: mode_file(:)
  type(String), allocatable :: line(:)
  integer                   :: no_atoms
  integer                   :: i
  
  mode_file = read_lines(filename)
  no_atoms = size(mode_file)-3
  call new(this,no_atoms)
  
  line = split(mode_file(1))
  this%frequency = dble(line(2))
  
  line = split(mode_file(2))
  this%soft_mode = line(3)=='T'
  
  do i=1,no_atoms
    line = split(mode_file(3+i))
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
  call print_line(mode_file, 'Frequency: '//this%frequency)
  call print_line(mode_file, 'Soft mode: '//this%soft_mode)
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
    
    ! Calculate q.R and exp(i q.R).
    qr = qpoint%qpoint * supercell%rvectors(supercell%atom_to_rvec(i))
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
end module
