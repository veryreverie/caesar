module mapping_module
  implicit none
  
  ! holds the contents of mapping.dat, and a few derived quantities
  type MappingData
    integer :: max
    integer :: first
    integer :: last
    integer :: count
    integer :: mid
  end type
  
  interface read_mapping_file
    module procedure read_mapping_file_character
    module procedure read_mapping_file_String
  end interface
  
contains

! reads a file ('mapping.dat'), and returns a MappingData
function read_mapping_file_character(filename) result(this)
  use file_module
  implicit none
  
  character(*), intent(in) :: filename
  type(MappingData)        :: this
  
  integer :: mapping_file
  
  mapping_file = open_read_file(filename)
  read(mapping_file,*) this%max
  read(mapping_file,*) this%first, this%last
  this%count = this%last-this%first+1
  this%mid = (this%count-1)/2+1
end function

function read_mapping_file_String(filename) result(this)
  use string_module
  implicit none
  
  type(String), intent(in) :: filename
  type(MappingData)        :: this
  
  this = read_mapping_file(char(filename))
end function
end module
