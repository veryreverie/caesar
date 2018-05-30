! ======================================================================
! The abstract potential type.
! ======================================================================
! See the example module in potential_pointer.f90 for how to use this type.
module potential_module
  use common_module
  
  use degeneracy_module
  use coupled_subspaces_module
  use degenerate_symmetry_module
  implicit none
  
  private
  
  public :: PotentialData
  public :: WriteLambda
  
  type, abstract, extends(NoDefaultConstructor) :: PotentialData
  contains
    procedure(generate_sampling_points_PotentialData), public, deferred :: &
       & generate_sampling_points
  end type
  
  interface
    subroutine WriteLambda(structure,directory)
      import StructureData
      import String
      implicit none
      
      type(StructureData), intent(in) :: structure
      type(String),        intent(in) :: directory
    end subroutine
  end interface
  
  abstract interface
    subroutine generate_sampling_points_PotentialData(this,                   &
       & sampling_points_dir,structure,symmetry_precision,complex_modes,      &
       & real_modes,qpoints,degenerate_subspaces,degenerate_symmetries,       &
       & coupled_subspaces,vscf_basis_functions_only,                         &
       & maximum_weighted_displacement,frequency_of_max_displacement,logfile, &
       & write_lambda)
      import dp
      import PotentialData
      import String
      import StructureData
      import CoupledSubspaces
      import RealMode
      import ComplexMode
      import QpointData
      import DegenerateModes
      import DegenerateSymmetry
      import CoupledSubspaces
      import OFile
      implicit none
      
      class(PotentialData),     intent(inout) :: this
      type(String),             intent(in)    :: sampling_points_dir
      type(StructureData),      intent(in)    :: structure
      real(dp),                 intent(in)    :: symmetry_precision
      type(ComplexMode),        intent(in)    :: complex_modes(:)
      type(RealMode),           intent(in)    :: real_modes(:)
      type(QpointData),         intent(in)    :: qpoints(:)
      type(DegenerateModes),    intent(in)    :: degenerate_subspaces(:)
      type(DegenerateSymmetry), intent(in)    :: degenerate_symmetries(:)
      type(CoupledSubspaces),   intent(in)    :: coupled_subspaces(:)
      logical,                  intent(in)    :: vscf_basis_functions_only
      real(dp),                 intent(in)    :: maximum_weighted_displacement
      real(dp),                 intent(in)    :: frequency_of_max_displacement
      type(OFile),              intent(inout) :: logfile
      procedure(WriteLambda)                  :: write_lambda
    end subroutine
  end interface
contains
end module
