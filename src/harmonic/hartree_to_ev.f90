! A Hartree to eV calculator
module caesar_hartree_to_ev_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_hartree_to_ev
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_hartree_to_ev() 
    end subroutine
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
