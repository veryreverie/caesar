module mapping_module
  implicit none
  
  ! holds the contents of mapping.dat, and a few derived quantities
  type Mapping
    integer :: max
    integer :: first
    integer :: last
    integer :: count
    integer :: mid
  end type
  
contains

! reads a file ('mapping.dat'), and returns a Mapping
function read_mapping(input) result(output)
  implicit none
  
  character(*), intent(in) :: input
  type(Mapping)            :: output
  
  open(100, file=input, status='old', action='read')
  read(100,*) output%max
  read(100,*) output%first, output%last
  close(100)
  output%count = output%last-output%first+1
  output%mid = (output%count-1)/2+1
end function

end module
