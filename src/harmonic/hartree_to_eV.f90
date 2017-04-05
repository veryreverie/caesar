! A Hartree to eV calculator
module hartree_to_eV_module
  implicit none
contains

subroutine hartree_to_eV()
  use constants, only : dp, ev_per_hartree
  use string_module
  use file_module
  implicit none
  
  real(dp) :: input
  
  call print_line('Input energy in Hartree:')
  input = dble(read_line_from_user())
  call print_line(input*ev_per_hartree)
end subroutine

end module
