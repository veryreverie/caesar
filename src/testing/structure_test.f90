! ======================================================================
! Provides helper functions for testing StructureData.
! ======================================================================
module caesar_structure_test_module
  use caesar_common_module
  implicit none
  
  interface
    ! ----------------------------------------------------------------------
    ! Writes out structure.dat in a manner compatible with old Caesar.
    ! ----------------------------------------------------------------------
    module subroutine write_old_structure_file(structure,filename) 
      type(StructureData), intent(in) :: structure
      type(String),        intent(in) :: filename
    end subroutine
  end interface
end module
