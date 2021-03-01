! ======================================================================
! Calculates observables under the VSCF approximation.
! ======================================================================
! Should be run after calculate_potential.
module caesar_calculate_anharmonic_observables_module
  use caesar_common_module 
  use caesar_states_module
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  
  use caesar_generate_subspace_potentials_module
  use caesar_vscf_module
  use caesar_interpolation_module
  implicit none
  
  private
  
  public :: startup_calculate_anharmonic_observables
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_calculate_anharmonic_observables() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! The main module subroutine.
    ! ----------------------------------------------------------------------
    module subroutine calculate_anharmonic_observables_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
