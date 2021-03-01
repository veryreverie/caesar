! ======================================================================
! Maps the potential along each normal mode.
! ======================================================================
module caesar_map_modes_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  
  use caesar_mode_map_module
  implicit none
  
  private
  
  public :: startup_map_modes
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_map_modes() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine map_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
