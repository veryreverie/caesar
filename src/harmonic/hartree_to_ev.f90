! A Hartree to eV calculator
module caesar_hartree_to_ev_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: hartree_to_ev_mode
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module function hartree_to_ev_mode() result(output)
      type(ProgramMode) :: output
    end function
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine hartree_to_ev_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
