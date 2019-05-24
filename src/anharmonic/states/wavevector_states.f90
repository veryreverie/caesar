! ======================================================================
! A set of states at a given wavevector.
! See wavevector_basis.f90 for more information.
! ======================================================================
module wavevector_states_module
  use common_module
  
  use wavevector_state_module
  implicit none
  
  private
  
  public :: WavevectorStates
  
  type, extends(NoDefaultConstructor) :: WavevectorStates
    type(WavevectorState), allocatable :: states(:)
    real(dp),              allocatable :: energies(:)
  end type
  
  interface WavevectorStates
    module procedure new_WavevectorStates
  end interface
contains

! Constructor.
function new_WavevectorStates(states,energies) result(this)
  implicit none
  
  type(WavevectorState), intent(in) :: states(:)
  real(dp),              intent(in) :: energies(:)
  type(WavevectorStates)            :: this
  
  this%states = states
  this%energies = energies
end function
end module
