! ======================================================================
! Forces in mass-weighted cartesian co-ordinates.
! ======================================================================
module mass_weighted_forces_submodule
  use utils_module
  implicit none
  
  private
  
  public :: MassWeightedForces
  
  type :: MassWeightedForces
    type(RealVector), allocatable :: forces(:)
  end type
contains
end module
