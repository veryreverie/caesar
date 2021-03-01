! ======================================================================
! Generate the set of q-points of the primitive structure which correspond
!    to G-vectors of a supercell.
! ======================================================================
module caesar_generate_qpoints_module
  use caesar_utils_module
  
  use caesar_structure_data_module
  use caesar_qpoint_module
  implicit none
  
  private
  
  public :: generate_qpoints
  
  interface
    module function generate_qpoints(large_supercell) result(output) 
      type(StructureData), intent(in) :: large_supercell
      type(QpointData), allocatable   :: output(:)
    end function
  end interface
end module
