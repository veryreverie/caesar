! ======================================================================
! The abstract potential type.
! ======================================================================
! See the example module in potential_pointer.f90 for how to use this type.
module potential_module
  use common_module
  
  use coupled_subspaces_module
  implicit none
  
  private
  
  public :: PotentialData
  
  type, abstract, extends(NoDefaultConstructor) :: PotentialData
  contains
    procedure(generate_sampling_points_PotentialData), public, deferred :: &
       & generate_sampling_points
  end type
  
  abstract interface
    subroutine generate_sampling_points_PotentialData(this,coupled_subspaces)
      import PotentialData
      import CoupledSubspaces
      implicit none
      
      class(PotentialData),   intent(inout) :: this
      type(CoupledSubspaces), intent(in)    :: coupled_subspaces(:)
    end subroutine
  end interface
contains
end module
