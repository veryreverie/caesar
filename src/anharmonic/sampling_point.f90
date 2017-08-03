! ======================================================================
! Holds information about points to be sampled.
! ======================================================================
module sampling_point_module
  use constants_module, only : dp
  use string_module
  use io_module
  
  type :: SamplingPoint
    integer, allocatable :: coordinates(:)
  end type
  
contains

subroutine write_sampling_points_file(this, filename)
  implicit none
  
  type(SamplingPoint), intent(in) :: this(:)
  type(String),        intent(in) :: filename
  
  integer :: sample_file
  integer :: i
  
  sample_file = open_write_file(filename)
  call print_line(sample_file,'! Sampling points in normal-mode co-ordinates.')
  call print_line(sample_file,'! Each co-ordinate is in multiples of delta.')
  do i=1,size(this)
    call print_line(sample_file,this(i)%coordinates)
  enddo
  close(sample_file)
end subroutine

function read_sampling_points_file(filename) result(this)
  implicit none
  
  type(String), intent(in)         :: filename
  type(SamplingPoint), allocatable :: this(:)
  
  type(String), allocatable :: sample_file(:)
  integer                   :: i,ialloc
  
  sample_file = read_lines(filename)
  
  allocate(this(size(sample_file)-2), stat=ialloc); call err(ialloc)
  
  do i=3,size(sample_file)
    this(i-2)%coordinates = int(split(sample_file(i)))
  enddo
end function
end module
