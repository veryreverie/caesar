! A Hartree to eV calculator
module hartree_to_eV_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function hartree_to_eV_keywords() result(keywords)
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(0)
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine hartree_to_eV(arguments)
  use constants_module, only : ev_per_hartree
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  real(dp) :: input
  
  call print_line('Input energy in Hartree:')
  input = dble(read_line_from_user())
  call print_line(input*ev_per_hartree)
end subroutine

end module
