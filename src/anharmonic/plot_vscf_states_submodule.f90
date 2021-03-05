submodule (caesar_plot_vscf_states_module) caesar_plot_vscf_states_submodule
  use caesar_anharmonic_module
contains

module procedure plot_vscf_states_mode
  output%mode_name = 'plot_vscf_states'
  output%description = 'plots the wavefunctions of the VSCF states, &
     &calculated by calculate_anharmonic_observables.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_vscf_states_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_vscf_states_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_states.py'),python_path)
end procedure
end submodule
