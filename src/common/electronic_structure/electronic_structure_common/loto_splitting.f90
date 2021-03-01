! ======================================================================
! Calculate the LO/TO correction to the energy and forces
!    at a given displacement,
!    using the electric field response calculated by DFPT.
! ======================================================================
module caesar_loto_splitting_module
  use caesar_utils_module
  
  use caesar_normal_mode_module
  use caesar_structure_module
  use caesar_electronic_structure_data_module
  implicit none
  
  private
  
  public :: LotoCorrection
  public :: calculate_loto_correction
  public :: dynamical_matrix_correction
  
  type, extends(NoDefaultConstructor) :: LotoCorrection
    real(dp)             :: energy
    type(CartesianForce) :: forces
  end type
  
  interface LotoCorrection
    ! Constructor.
    module function new_LotoCorrection(energy,forces) result(this) 
      real(dp),             intent(in) :: energy
      type(CartesianForce), intent(in) :: forces
      type(LotoCorrection)             :: this
    end function
  end interface
  
  interface LotoCorrection
    ! N.B. The force correction is only valid under the assumption that
    !    the Born effective charges and the permitivity are constants.
    module function new_LotoCorrection_LinearResponse(linear_response, &
       & loto_direction,displacement,structure) result(this) 
      type(LinearResponse),        intent(in) :: linear_response
      type(FractionVector),        intent(in) :: loto_direction
      type(CartesianDisplacement), intent(in) :: displacement
      type(StructureData),         intent(in) :: structure
      type(LotoCorrection)                    :: this
    end function
  end interface
  
  interface
    ! Calculate the dynamical matrix correction.
    module function dynamical_matrix_correction(linear_response, &
       & loto_direction,structure) result(output) 
      type(LinearResponse), intent(in) :: linear_response
      type(FractionVector), intent(in) :: loto_direction
      type(StructureData),  intent(in) :: structure
      type(RealMatrix), allocatable    :: output(:,:)
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Uses calculated LO/TO correction to update electronic structure output.
    ! ----------------------------------------------------------------------
    module function calculate_loto_correction(electronic_structure, &
       & loto_correction) result(output) 
      type(ElectronicStructure), intent(in) :: electronic_structure
      type(LotoCorrection),      intent(in) :: loto_correction
      type(ElectronicStructure)             :: output
    end function
  end interface
end module
