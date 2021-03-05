! ======================================================================
! Maps out the potential at finite displacements along multiple modes.
! ======================================================================
module caesar_map_potential_module
  use caesar_common_module
  use caesar_harmonic_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  implicit none
  
  private
  
  public :: map_potential_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function map_potential_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine map_potential_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
