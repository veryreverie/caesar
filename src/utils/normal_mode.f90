module normal_mode_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  type NormalMode
    logical               :: soft_mode
    real(dp)              :: frequency
    real(dp), allocatable :: displacements(:,:)
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
  
  allocate(this%displacements(3,no_atoms), stat=ialloc); call err(ialloc)
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
    this%displacements(:,i) = dble(line)
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
  do i=1,size(this%displacements,2)
    call print_line(mode_file, this%displacements(:,i))
  enddo
  close(mode_file)
end subroutine
end module
