! A Hartree to eV calculator
module hartree_to_ev_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function hartree_to_ev_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'energy_in_hartree',                                        &
  &               'energy_in_hartree is the input energy, in Hartree atomic &
  &units.') ]
end function

function hartree_to_ev_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'hartree_to_ev'
  output%description = 'Converts energy in eV to energy in Hartree.'
  output%keywords = hartree_to_ev_keywords()
  output%main_subroutine => hartree_to_ev
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine hartree_to_ev(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  real(dp) :: input
  
  input = dble(arguments%value('energy_in_hartree'))
  call print_line(input*ev_per_hartree)
end subroutine

end module
