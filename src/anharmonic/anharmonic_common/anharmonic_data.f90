! ======================================================================
! A storage class, for all the variables which are passed to the
!    various potential methods.
! ======================================================================
module anharmonic_data_module
  use common_module
  
  use degenerate_symmetry_module
  use subspace_coupling_module
  implicit none
  
  private
  
  public :: AnharmonicData
  
  type, extends(NoDefaultConstructor) :: AnharmonicData
    type(StructureData)                   :: structure
    type(StructureData)                   :: anharmonic_supercell
    type(QpointData),         allocatable :: qpoints(:)
    type(ComplexMode),        allocatable :: complex_modes(:)
    type(RealMode),           allocatable :: real_modes(:)
    type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
    type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
    type(SubspaceCoupling),   allocatable :: subspace_couplings(:)
    logical                               :: vscf_basis_functions_only
    real(dp)                              :: maximum_weighted_displacement
    real(dp)                              :: frequency_of_max_displacement
  end type
  
  interface AnharmonicData
    module procedure new_AnharmonicData
  end interface
contains

function new_AnharmonicData(structure,anharmonic_supercell,qpoints,       &
   & complex_modes,real_modes,degenerate_subspaces,degenerate_symmetries, &
   & subspace_couplings,vscf_basis_functions_only,                        &
   & maximum_weighted_displacement,frequency_of_max_displacement) result(this)
  implicit none
  
  type(StructureData),      intent(in) :: structure
  type(StructureData),      intent(in) :: anharmonic_supercell
  type(QpointData),         intent(in) :: qpoints(:)
  type(ComplexMode),        intent(in) :: complex_modes(:)
  type(RealMode),           intent(in) :: real_modes(:)
  type(DegenerateSubspace), intent(in) :: degenerate_subspaces(:)
  type(DegenerateSymmetry), intent(in) :: degenerate_symmetries(:)
  type(SubspaceCoupling),   intent(in) :: subspace_couplings(:)
  logical,                  intent(in) :: vscf_basis_functions_only
  real(dp),                 intent(in) :: maximum_weighted_displacement
  real(dp),                 intent(in) :: frequency_of_max_displacement
  type(AnharmonicData)                 :: this
  
  this%structure                     = structure
  this%anharmonic_supercell          = anharmonic_supercell
  this%qpoints                       = qpoints
  this%complex_modes                 = complex_modes
  this%real_modes                    = real_modes
  this%degenerate_subspaces          = degenerate_subspaces
  this%degenerate_symmetries         = degenerate_symmetries
  this%subspace_couplings            = subspace_couplings
  this%vscf_basis_functions_only     = vscf_basis_functions_only
  this%maximum_weighted_displacement = maximum_weighted_displacement
  this%frequency_of_max_displacement = frequency_of_max_displacement
end function
end module
