submodule (caesar_plot_cutoff_and_kpoints_module) caesar_plot_cutoff_and_kpoints_submodule
  use caesar_dft_module
contains

module procedure plot_cutoff_and_kpoints
  output%mode_name = 'plot_cutoff_and_kpoints'
  output%description = 'Plots the output of coverge_cutoff_and_kpoints'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_cutoff_and_kpoints_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_cutoff_and_kpoints_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_cutoff_and_kpoints.py'),python_path)
end procedure
end submodule
