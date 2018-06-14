! ======================================================================
! Plots the phonon density of states and dispersion calculated by
!    calculate_harmonic_observables.
! ======================================================================
module plot_dos_and_dispersion_module
  use common_module
  implicit none
  
  private
  
  public :: plot_dos_and_dispersion
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_dos_and_dispersion() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_dos_and_dispersion'
  output%description = 'Plots the phonon density of states and dispersion &
     &calculated by calculate_harmonic_observables.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_dos_and_dispersion_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_dos_and_dispersion_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_dos_and_dispersion.py'),python_path)
end subroutine
end module
