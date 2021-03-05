! ======================================================================
! Calculates observables under the VSCF approximation.
! ======================================================================
! Should be run after calculate_potential.
module caesar_calculate_anharmonic_observables_module
  use caesar_common_module 
  use caesar_harmonic_module
  
  use caesar_states_module
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  
  use caesar_generate_subspace_potentials_module
  use caesar_vscf_module
  use caesar_interpolation_module
  implicit none
  
  private
  
  public :: calculate_anharmonic_observables_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function calculate_anharmonic_observables_mode() result(output)
      type(ProgramMode) :: output
    end function
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
