submodule (caesar_plot_vscf_convergence_module) caesar_plot_vscf_convergence_submodule
  use caesar_anharmonic_module
contains

module procedure startup_plot_vscf_convergence
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_vscf_convergence'
  mode%description = 'Plots the convergence of the VSCF scheme. This mode &
     &should be called from anharmonic_observables/*/temperature_*.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_vscf_convergence_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end procedure

module procedure plot_vscf_convergence_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_convergence.py'),python_path)
end procedure
end submodule
