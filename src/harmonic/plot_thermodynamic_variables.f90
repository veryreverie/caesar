! ======================================================================
! Plots the thermodynamic variables calculated by
!    calculate_harmonic_observables.
! ======================================================================
module plot_thermodynamic_variables_module
  use common_module
  implicit none
  
  private
  
  public :: plot_thermodynamic_variables
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_thermodynamic_variables() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_thermodynamic_variables'
  output%description = 'Plots the thermodynamic variables &
     &calculated by calculate_harmonic_observables.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_thermodynamic_variables_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_thermodynamic_variables_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_thermodynamic_variables.py'),python_path)
end subroutine
end module
