submodule (caesar_plot_potential_map_module) caesar_plot_potential_map_submodule
  use caesar_anharmonic_module
contains

module procedure plot_potential_map_mode
  output%mode_name = 'plot_potential_map'
  output%description = 'Plots the mapping of the anharmonic potential &
     &calculated by map_potential.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_potential_map_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_potential_map_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_potential_map.py'),python_path)
end procedure
end submodule
