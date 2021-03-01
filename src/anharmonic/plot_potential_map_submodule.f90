submodule (caesar_plot_potential_map_module) caesar_plot_potential_map_submodule
  use caesar_anharmonic_module
contains

module procedure startup_plot_potential_map
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_potential_map'
  mode%description = 'Plots the mapping of the anharmonic potential &
     &calculated by map_potential.'
  mode%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_potential_map_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end procedure

module procedure plot_potential_map_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_potential_map.py'),python_path)
end procedure
end submodule
