! ======================================================================
! Reads and writes structure.dat files from StructureData.
! ======================================================================
module structure_file_module
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: read_structure_file
  public :: write_structure_file
contains

function read_structure_file(filename) result(this)
  implicit none
  
  type(String), intent(in) :: filename
  type(StructureData)      :: this
  
  type(IFile) :: structure_file
  
  structure_file = IFile(filename)
  this = StructureData(structure_file%lines())
end function

subroutine write_structure_file(this,filename)
  implicit none
  
  type(StructureData), intent(in) :: this
  type(String),        intent(in) :: filename
  
  type(OFile) :: structure_file
  
  structure_file = OFile(filename)
  call structure_file%print_lines(this)
end subroutine
end module
