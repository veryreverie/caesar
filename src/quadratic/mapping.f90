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
function read_mapping(mapping_file_unit) result(output)
  implicit none
  
  integer, intent(in) :: mapping_file_unit
  type(Mapping)       :: output
  
  read(mapping_file_unit,*) output%max
  read(mapping_file_unit,*) output%first, output%last
  output%count = output%last-output%first+1
  output%mid = (output%count-1)/2+1
end function

end module
