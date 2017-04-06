! A Hartree to eV calculator
module hartree_to_eV_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

subroutine hartree_to_eV()
  use constants_module, only : ev_per_hartree
  implicit none
  
  real(dp) :: input
  
  call print_line('Input energy in Hartree:')
  input = dble(read_line_from_user())
  call print_line(input*ev_per_hartree)
end subroutine

end module
