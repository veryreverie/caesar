submodule (caesar_plot_thermodynamic_variables_module) caesar_plot_thermodynamic_variables_submodule
  use caesar_harmonic_module
contains

module procedure plot_thermodynamic_variables_mode
  output%mode_name = 'plot_thermodynamic_variables'
  output%description = 'Plots the thermodynamic variables &
     &calculated by calculate_harmonic_observables or &
     & calculate_anharmonic_observables. Should be called from within the &
     &harmonic_observables or anharmonic_observables directory. The -d flag &
     &may be useful for this.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_thermodynamic_variables_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_thermodynamic_variables_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_thermodynamic_variables.py'),python_path)
end procedure
end submodule
