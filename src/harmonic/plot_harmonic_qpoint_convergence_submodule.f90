submodule (caesar_plot_harmonic_qpoint_convergence_module) caesar_plot_harmonic_qpoint_convergence_submodule
  use caesar_harmonic_module
contains

module procedure startup_plot_harmonic_qpoint_convergence
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_harmonic_qpoint_convergence'
  mode%description = 'Plots the convergence of the harmonic free energy &
     &calculation w/r/t the q-point grid, as calculated by &
     &converge_harmonic_qpoints.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_harmonic_qpoint_convergence_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end procedure

module procedure plot_harmonic_qpoint_convergence_subroutine
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_harmonic_qpoint_convergence.py'),python_path)
end procedure
end submodule
