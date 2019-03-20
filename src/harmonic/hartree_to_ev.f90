! A Hartree to eV calculator
module hartree_to_ev_module
  use common_module
  implicit none
  
  private
  
  public :: startup_hartree_to_ev
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_hartree_to_ev()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'hartree_to_ev'
  mode%description = 'Converts energy in eV to energy in Hartree.'
  mode%keywords = [                                                        &
  & KeywordData( 'energy_in_hartree',                                      &
  &              'energy_in_hartree is the input energy, in Hartree atomic &
  &units.') ]
  mode%main_subroutine => hartree_to_ev_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine hartree_to_ev_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  real(dp) :: input
  
  input = dble(arguments%value('energy_in_hartree'))
  call print_line(input*EV_PER_HARTREE)
end subroutine

end module
