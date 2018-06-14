! A Hartree to eV calculator
module hartree_to_ev_module
  use common_module
  implicit none
  
  private
  
  public :: hartree_to_ev
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function hartree_to_ev() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'hartree_to_ev'
  output%description = 'Converts energy in eV to energy in Hartree.'
  output%keywords = [                                                      &
  & KeywordData( 'energy_in_hartree',                                      &
  &              'energy_in_hartree is the input energy, in Hartree atomic &
  &units.') ]
  output%main_subroutine => hartree_to_ev_subroutine
end function

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
