! ======================================================================
! Plots the VSCF potential along each mode.
! ======================================================================
module plot_vscf_states_module
  use common_module
  implicit none
  
  private
  
  public :: plot_vscf_states
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_vscf_states() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_vscf_states'
  output%description = 'plots the wavefunctions of the VSCF states, &
     &calculated by calculate_anharmonic_observables.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_vscf_states_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_vscf_states_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_states.py'),python_path)
end subroutine
end module
