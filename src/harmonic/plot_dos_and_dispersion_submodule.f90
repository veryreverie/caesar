submodule (caesar_plot_dos_and_dispersion_module) caesar_plot_dos_and_dispersion_submodule
  use caesar_harmonic_module
contains

module procedure plot_dos_and_dispersion_mode
  output%mode_name = 'plot_dos_and_dispersion'
  output%description = 'Plots the phonon density of states and dispersion &
     &calculated by calculate_harmonic_observables or &
     &calculate_anharmonic_observables. Should be run from within the &
     &harmonic_observables, anharmonic_observables/* or &
     &anharmonic_observables/*/temperature_* directory. The -d flag may be &
     &useful for this.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_dos_and_dispersion_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_dos_and_dispersion_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_dos_and_dispersion.py'),python_path)
end procedure
end submodule
