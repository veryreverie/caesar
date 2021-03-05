! ======================================================================
! Maps the VSCF potential along each normal mode.
! ======================================================================
module caesar_map_vscf_modes_module
  use caesar_common_module
  
  use caesar_states_module
  use caesar_anharmonic_common_module
  use caesar_potentials_module
  
  use caesar_mode_map_module
  implicit none
  
  private
  
  public :: map_vscf_modes_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function map_vscf_modes_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine map_vscf_modes_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
