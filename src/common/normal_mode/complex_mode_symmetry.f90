! ======================================================================
! Methods for converting symmetries into complex normal mode co-ordinates.
! ======================================================================
module caesar_complex_mode_symmetry_module
  use caesar_utils_module
  
  use caesar_structure_module
  
  use caesar_complex_mode_module
  implicit none
  
  private
  
  public :: calculate_symmetry_in_normal_coordinates
  
  interface
    ! ----------------------------------------------------------------------
    ! Calculates a symmetry in normal mode co-ordinates.
    ! ----------------------------------------------------------------------
    module function calculate_symmetry_in_normal_coordinates(modes,qpoints, &
       & symmetry) result(output) 
      type(ComplexMode),      intent(in)    :: modes(:)
      type(QpointData),       intent(in)    :: qpoints(:)
      type(SymmetryOperator), intent(in)    :: symmetry
      type(ComplexMatrix)                   :: output
    end function
  end interface
end module
