! ======================================================================
! Plots the convergence of the harmonic free energy calculation
!    w/r/t the q-point grid, as calculated by converge_harmonic_qpoints.
! ======================================================================
module caesar_plot_harmonic_qpoint_convergence_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_harmonic_qpoint_convergence
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_harmonic_qpoint_convergence()
  implicit none
  
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
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_harmonic_qpoint_convergence_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_harmonic_qpoint_convergence.py'),python_path)
end subroutine
end module
