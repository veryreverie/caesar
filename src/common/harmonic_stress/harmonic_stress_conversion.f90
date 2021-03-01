! ======================================================================
! Conversions between cartesian and normal-mode representations of
!    harmonic stress.
! ======================================================================
module caesar_harmonic_stress_conversion_module
  use caesar_utils_module
  use caesar_structure_module
  use caesar_normal_mode_module
  use caesar_dynamical_matrices_module
  
  use caesar_stress_hessian_module
  use caesar_stress_dynamical_matrix_module
  implicit none
  
  private
  
  public :: StressDynamicalMatrix
  public :: reconstruct_stress_hessian
  
  interface StressDynamicalMatrix
    ! Construct a stress dynamical matrix at an arbitrary q-point.
    module function new_StressDynamicalMatrix_interpolated(q,supercell, &
       & hessian,min_images) result(output) 
      type(RealVector),    intent(in)           :: q
      type(StructureData), intent(in)           :: supercell
      type(StressHessian), intent(in)           :: hessian
      type(MinImages),     intent(in), optional :: min_images(:,:)
      type(StressDynamicalMatrix)               :: output
    end function
  end interface
  
  interface
    module function reconstruct_stress_hessian(large_supercell,qpoints, &
       & dynamical_matrices,logfile) result(output) 
      type(StructureData),         intent(in)    :: large_supercell
      type(QpointData),            intent(in)    :: qpoints(:)
      type(StressDynamicalMatrix), intent(in)    :: dynamical_matrices(:)
      type(OFile),                 intent(inout) :: logfile
      type(StressHessian)                        :: output
    end function
  end interface
end module
