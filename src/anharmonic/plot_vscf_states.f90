! ======================================================================
! Plots the VSCF potential along each mode.
! ======================================================================
module caesar_plot_vscf_states_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_vscf_states
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_vscf_states()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_vscf_states'
  mode%description = 'plots the wavefunctions of the VSCF states, &
     &calculated by calculate_anharmonic_observables.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_vscf_states_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

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
