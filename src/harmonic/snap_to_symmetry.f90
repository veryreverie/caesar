! ======================================================================
! Uses spglib to snap a structure to a symmetry.
! ======================================================================
module caesar_snap_to_symmetry_module
  use caesar_common_module
  
  use caesar_generate_supercells_module
  implicit none
  
  private
  
  public :: startup_snap_to_symmetry
  
  interface
    ! ----------------------------------------------------------------------
    ! Generate keywords and helptext.
    ! ----------------------------------------------------------------------
    module subroutine startup_snap_to_symmetry() 
    end subroutine
  end interface
  
  interface
    ! ----------------------------------------------------------------------
    ! Main program.
    ! ----------------------------------------------------------------------
    module subroutine snap_to_symmetry_subroutine(arguments) 
      type(Dictionary), intent(in) :: arguments
    end subroutine
  end interface
end module
