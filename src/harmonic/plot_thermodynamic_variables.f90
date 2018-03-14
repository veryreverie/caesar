! ======================================================================
! Plots the thermodynamic variables calculated by
!    calculate_harmonic_observables.
! ======================================================================
module plot_thermodynamic_variables_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_thermodynamic_variables_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                             &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
end function

function plot_thermodynamic_variables_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_thermodynamic_variables'
  output%description = 'Plots the thermodynamic variables &
     &calculated by calculate_harmonic_observables.'
  output%keywords = plot_thermodynamic_variables_keywords()
  output%main_subroutine => plot_thermodynamic_variables
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_thermodynamic_variables(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_thermodynamic_variables.py'),python_path)
end subroutine
end module
