! ======================================================================
! Interpolate an anharmonic dynamical matrix to the same q-point density
!    as a given harmonic dynamical matrix.
! ======================================================================
module caesar_interpolation_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  
  private
  
  public :: interpolate_dynamical_matrices
  public :: calculate_difference_dynamical_matrices
  public :: calculate_interpolated_hessian
  public :: calculate_interpolated_stress_hessian
  
  interface
    ! Takes a set of dynamical matrices defined on a coarse q-point grid,
    !    and interpolates to give the dynamical matrices on a fine q-point grid.
    module function interpolate_dynamical_matrices(coarse_dynamical_matrices,coarse_qpoints,coarse_supercell,coarse_min_images,fine_qpoints) result(output) 
      type(DynamicalMatrix), intent(in)  :: coarse_dynamical_matrices(:)
      type(QpointData),      intent(in)  :: coarse_qpoints(:)
      type(StructureData),   intent(in)  :: coarse_supercell
      type(MinImages),       intent(in)  :: coarse_min_images(:,:)
      type(QpointData),      intent(in)  :: fine_qpoints(:)
      type(DynamicalMatrix), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! When dynamical matrices are interpolated from a fine q-point grid
    !    onto a coarse q-point grid, some information may be lost.
    ! This module function takes two inputs:
    !    - a set of dynamical matrices defined on a fine q-point grid.
    !    - a coarse q-point grid.
    ! The module function separates the dynamical matrices into the sum of two parts: &
    !    - the part which is preserved when the dynamical matrices are interpolated   &
    !         onto the coarse q-point grid and back again.
    !    - the part which is lost when the dynamical matrices are interpolated
    !         onto the coarse q-point grid and back again.
    ! The module function then returns only the contribution to the dynamical matrices
    !    which is lost under interpolation.
    module function calculate_difference_dynamical_matrices(fine_dynamical_matrices,fine_qpoints,fine_supercell,fine_min_images,coarse_qpoints,coarse_supercell,coarse_min_images) result(output) 
      type(DynamicalMatrix), intent(in)  :: fine_dynamical_matrices(:)
      type(QpointData),      intent(in)  :: fine_qpoints(:)
      type(StructureData),   intent(in)  :: fine_supercell
      type(MinImages),       intent(in)  :: fine_min_images(:,:)
      type(QpointData),      intent(in)  :: coarse_qpoints(:)
      type(StructureData),   intent(in)  :: coarse_supercell
      type(MinImages),       intent(in)  :: coarse_min_images(:,:)
      type(DynamicalMatrix), allocatable :: output(:)
    end function
  end interface
  
  interface
    ! Takes anharmonic dynamical matrices on a coarse q-point grid,
    !    and harmonic dynamical matrices on a fine grid,
    !    and returns a fine-grid Hessian containing both.
    module function calculate_interpolated_hessian(anharmonic_supercell, &
       & anharmonic_qpoints,anharmonic_dynamical_matrices,               &
       & anharmonic_min_images,harmonic_supercell,harmonic_qpoints,      &
       & harmonic_dynamical_matrices,harmonic_min_images) result(output) 
      type(StructureData),   intent(in) :: anharmonic_supercell
      type(QpointData),      intent(in) :: anharmonic_qpoints(:)
      type(DynamicalMatrix), intent(in) :: anharmonic_dynamical_matrices(:)
      type(MinImages),       intent(in) :: anharmonic_min_images(:,:)
      type(StructureData),   intent(in) :: harmonic_supercell
      type(QpointData),      intent(in) :: harmonic_qpoints(:)
      type(DynamicalMatrix), intent(in) :: harmonic_dynamical_matrices(:)
      type(MinImages),       intent(in) :: harmonic_min_images(:,:)
      type(CartesianHessian)            :: output
    end function
  end interface
  
  interface
    module function calculate_interpolated_stress_hessian(anharmonic_supercell,anharmonic_qpoints,anharmonic_dynamical_matrices,anharmonic_min_images,harmonic_supercell,harmonic_qpoints,harmonic_dynamical_matrices,harmonic_min_images) result(output) 
      type(StructureData),         intent(in) :: anharmonic_supercell
      type(QpointData),            intent(in) :: anharmonic_qpoints(:)
      type(StressDynamicalMatrix), intent(in) :: anharmonic_dynamical_matrices(:)
      type(MinImages),             intent(in) :: anharmonic_min_images(:,:)
      type(StructureData),         intent(in) :: harmonic_supercell
      type(QpointData),            intent(in) :: harmonic_qpoints(:)
      type(StressDynamicalMatrix), intent(in) :: harmonic_dynamical_matrices(:)
      type(MinImages),             intent(in) :: harmonic_min_images(:,:)
      type(StressHessian)                     :: output
    end function
  end interface
end module
