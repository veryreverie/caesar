submodule (caesar_plot_vscf_convergence_module) caesar_plot_vscf_convergence_submodule
  use caesar_anharmonic_module
contains

module procedure plot_vscf_convergence_mode
  output%mode_name = 'plot_vscf_convergence'
  output%description = 'Plots the convergence of the VSCF scheme. This mode &
     &should be called from anharmonic_observables/*/temperature_*.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_vscf_convergence_subroutine
  output%suppress_settings_file = .true.
end procedure

module procedure plot_vscf_convergence_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_convergence.py'),python_path)
end procedure
end submodule
