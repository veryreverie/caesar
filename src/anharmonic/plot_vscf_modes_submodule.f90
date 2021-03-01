submodule (caesar_plot_vscf_modes_module) caesar_plot_vscf_modes_submodule
  use caesar_anharmonic_module
contains

module procedure startup_plot_vscf_modes
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_vscf_modes'
  mode%description = 'Plots the mapping of the vscf modes produced by &
     &map_vscf_modes.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_vscf_modes_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end procedure

module procedure plot_vscf_modes_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_modes.py'),python_path)
end procedure
end submodule
