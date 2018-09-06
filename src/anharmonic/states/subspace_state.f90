! ======================================================================
! A state spanning a given subspace.
! ======================================================================
module subspace_state_module
  use common_module
  implicit none
  
  private
  
  public :: SubspaceState
  
  type, abstract, extends(Stringsable) :: SubspaceState
    integer :: subspace_id
  end type
end module
