submodule (caesar_plot_modes_module) caesar_plot_modes_submodule
  use caesar_anharmonic_module
contains

module procedure plot_modes_mode
  output%mode_name = 'plot_modes'
  output%description = 'Plots the modes mapped by map_modes.'
  output%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_modes_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_modes_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_modes.py'),python_path)
end procedure
end submodule
