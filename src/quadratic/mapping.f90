module mapping_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
  ! holds the contents of mapping.dat, and a few derived quantities
  type MappingData
    real(dp) :: max
    integer  :: first
    integer  :: last
    integer  :: count
    integer  :: mid
  end type
  
  interface read_mapping_file
    module procedure read_mapping_file_character
    module procedure read_mapping_file_String
  end interface
  
contains

! reads a file ('mapping.dat'), and returns a MappingData
function read_mapping_file_character(filename) result(this)
  implicit none
  
  character(*), intent(in) :: filename
  type(MappingData)        :: this
  
  type(String), allocatable :: mapping_file(:)
  type(String), allocatable :: line(:)
  
  mapping_file = read_lines(filename)
  this%max = dble(mapping_file(1))
  line = split(mapping_file(2))
  this%first = int(line(1))
  this%last = int(line(2))
  this%count = this%last-this%first+1
  this%mid = (this%count-1)/2+1
end function

function read_mapping_file_String(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(MappingData)        :: this
  
  this = read_mapping_file(char(filename))
end function
end module
