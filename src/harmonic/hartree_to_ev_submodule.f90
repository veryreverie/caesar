submodule (caesar_hartree_to_ev_module) caesar_hartree_to_ev_submodule
  use caesar_harmonic_module
contains

module procedure hartree_to_ev_mode
  output%mode_name = 'hartree_to_ev'
  output%description = 'Converts energy in eV to energy in Hartree.'
  output%keywords = [                                                        &
  & KeywordData( 'energy_in_hartree',                                      &
  &              'energy_in_hartree is the input energy, in Hartree atomic &
  &units.') ]
  output%main_subroutine => hartree_to_ev_subroutine
end procedure

module procedure hartree_to_ev_subroutine
  real(dp) :: input
  
  input = dble(arguments%value('energy_in_hartree'))
  call print_line(input*EV_PER_HARTREE)
end procedure
end submodule
