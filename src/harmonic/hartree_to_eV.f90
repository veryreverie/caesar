! A Hartree to eV calculator
module hartree_to_eV_module
  implicit none
contains

subroutine hartree_to_eV()
  use constants, only : dp, eV
  implicit none
  
  real(dp) :: input
  
  write(*,"(a)") "Input energy in Hartree:"
  read(*,*) input
  write(*,*) input*eV
end subroutine

end module
