! ======================================================================
! Updates basis function files from Dec2018 format to Jan2019 format.
! ======================================================================
module caesar_update_basis_functions_module
  use caesar_common_module
  
  use caesar_anharmonic_common_module
  implicit none
  
  private
  
  public :: startup_update_basis_functions
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_update_basis_functions() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine update_basis_functions_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
