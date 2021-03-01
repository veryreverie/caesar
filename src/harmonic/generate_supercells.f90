! ======================================================================
! Generates the supercells needed to simulate harmonic phonons at all q-points.
! ======================================================================
module caesar_generate_supercells_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: generate_supercells
  
  interface
    ! ----------------------------------------------------------------------
    ! Generates supercells.
    ! ----------------------------------------------------------------------
    module function generate_supercells(structure,qpoints, &
       & symmetry_precision,loto_direction) result(output) 
      ! Inputs.
      type(StructureData),  intent(in)           :: structure
      type(QpointData),     intent(in)           :: qpoints(:)
      real(dp),             intent(in)           :: symmetry_precision
      type(FractionVector), intent(in), optional :: loto_direction
      type(StructureData), allocatable           :: output(:)
    end function
  end interface
end module
