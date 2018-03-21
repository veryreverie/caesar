! ======================================================================
! A class for holding calculated energy and forces.
! ======================================================================
module electronic_structure_data_submodule
  use utils_module
  
  use structure_module
  implicit none
  
  private
  
  public :: ElectronicStructure
  
  type, extends(NoDefaultConstructor) :: ElectronicStructure
    real(dp)                      :: energy
    type(RealVector), allocatable :: forces(:)
  end type
  
  interface ElectronicStructure
    module procedure new_ElectronicStructure
  end interface
contains

function new_ElectronicStructure(energy,forces) result(this)
  implicit none
  
  real(dp),         intent(in) :: energy
  type(RealVector), intent(in) :: forces(:)
  type(ElectronicStructure)    :: this
  
  this%energy  = energy
  this%forces  = forces
end function
end module
