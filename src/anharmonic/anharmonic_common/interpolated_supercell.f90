! ======================================================================
! Generates a supercell with a given q-point grid,
!    along with the corresponding q-points and normal modes.
! ======================================================================
module caesar_interpolated_supercell_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: InterpolatedSupercell
  
  type, extends(NoDefaultConstructor) :: InterpolatedSupercell
    type(StructureData)                :: supercell
    type(QpointData),      allocatable :: qpoints(:)
    type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
    type(ComplexMode),     allocatable :: complex_modes(:)
    type(RealMode),        allocatable :: real_modes(:)
  end type
  
  interface InterpolatedSupercell
    module function new_InterpolatedSupercell(supercell,qpoints, &
       & dynamical_matrices,complex_modes,real_modes) result(this) 
      type(StructureData),   intent(in) :: supercell
      type(QpointData),      intent(in) :: qpoints(:)
      type(DynamicalMatrix), intent(in) :: dynamical_matrices(:)
      type(ComplexMode),     intent(in) :: complex_modes(:)
      type(RealMode),        intent(in) :: real_modes(:)
      type(InterpolatedSupercell)       :: this
    end function
  
    module function new_InterpolatedSupercell_interpolated(qpoint_grid, &
       & structure,harmonic_supercell,harmonic_qpoints,                 &
       & harmonic_dynamical_matrices,harmonic_complex_modes,logfile)    &
       & result(this) 
      integer,               intent(in)    :: qpoint_grid(:)
      type(StructureData),   intent(in)    :: structure
      type(StructureData),   intent(in)    :: harmonic_supercell
      type(QpointData),      intent(in)    :: harmonic_qpoints(:)
      type(DynamicalMatrix), intent(in)    :: harmonic_dynamical_matrices(:)
      type(ComplexMode),     intent(in)    :: harmonic_complex_modes(:,:)
      type(OFile),           intent(inout) :: logfile
      type(InterpolatedSupercell)          :: this
    end function
  end interface
end module
