submodule (caesar_plot_normal_modes_module) caesar_plot_normal_modes_submodule
  use caesar_harmonic_module
contains

module procedure plot_normal_modes_mode
  output%mode_name = 'plot_normal_modes'
  output%description = 'Plots the output of calculate_normal_modes. Should be &
     &run from within a qpoint_ directory. The -d flag may be useful for this.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_normal_modes_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_normal_modes_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_normal_modes.py'),python_path)
end procedure
end submodule
