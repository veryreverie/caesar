! ======================================================================
! Reads and writes structure.dat files from StructureData.
! ======================================================================
module caesar_structure_file_module
  use caesar_utils_module
  
  use caesar_structure_module
  implicit none
  
  private
  
  public :: read_structure_file
  public :: write_structure_file
  
  interface
    module function read_structure_file(filename) result(this) 
      type(String), intent(in) :: filename
      type(StructureData)      :: this
    end function
  end interface
  
  interface
    module subroutine write_structure_file(this,filename) 
      type(StructureData), intent(in) :: this
      type(String),        intent(in) :: filename
    end subroutine
  end interface
end module
