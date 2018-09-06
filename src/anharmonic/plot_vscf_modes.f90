! ======================================================================
! Plots the VSCF potential along each mode.
! ======================================================================
module plot_vscf_modes_module
  use common_module
  implicit none
  
  private
  
  public :: plot_vscf_modes
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_vscf_modes() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_vscf_modes'
  output%description = 'Plots the phonon density of states and dispersion &
     &calculated by calculate_harmonic_observables. Should be run from within &
     &a qpoint_ directory. The -d flag may be useful for this.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_vscf_modes_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_vscf_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_vscf_modes.py'),python_path)
end subroutine
end module
