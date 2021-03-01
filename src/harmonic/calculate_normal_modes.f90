! ======================================================================
! The third stage of Caesar.
! Uses the forces calculated previously to calculate harmonic normal modes.
! ======================================================================
! Calculates the set of phonon normal modes at each supercell G-vector.
! Calculates the real part of the non-mass-reduced polarisation vector, which
!    is the pattern of displacement corresponding to the normal mode.
module caesar_calculate_normal_modes_module
  use caesar_common_module
  
  use caesar_construct_supercell_hessian_module
  implicit none
  
  private
  
  public :: startup_calculate_normal_modes
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_calculate_normal_modes() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine calculate_normal_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
