module supercells_module
contains

function read_supercells(filename) result(this)
  use string_module
  use file_module
  implicit none
  
  type(String), intent(in) :: filename
  integer, allocatable     :: this(:,:,:)
  
  integer                   :: i,j
  type(String), allocatable :: file_contents(:)
  
  file_contents = read_lines(filename)
  
  allocate(this(3,3,size(file_contents)/4))
  
  do i=1,size(file_contents)/4
    do j=1,3
      this(j,:,i) = int(split(file_contents(4*(i-1)+j)))
    enddo
  enddo
end function
end module
