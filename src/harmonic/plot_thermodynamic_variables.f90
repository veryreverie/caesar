! ======================================================================
! Plots the thermodynamic variables calculated by
!    calculate_harmonic_observables or calculate_anharmonic_observables.
! ======================================================================
module caesar_plot_thermodynamic_variables_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_thermodynamic_variables
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_thermodynamic_variables()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_thermodynamic_variables'
  mode%description = 'Plots the thermodynamic variables &
     &calculated by calculate_harmonic_observables or &
     & calculate_anharmonic_observables. Should be called from within the &
     &harmonic_observables or anharmonic_observables directory. The -d flag &
     &may be useful for this.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_thermodynamic_variables_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

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
