! ======================================================================
! Plots the phonon density of states and dispersion calculated by
!    calculate_harmonic_observables or calculate_anharmonic_observables.
! ======================================================================
module plot_dos_and_dispersion_module
  use common_module
  implicit none
  
  private
  
  public :: startup_plot_dos_and_dispersion
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_dos_and_dispersion()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_dos_and_dispersion'
  mode%description = 'Plots the phonon density of states and dispersion &
     &calculated by calculate_harmonic_observables or &
     &calculate_anharmonic_observables. Should be run from within the &
     &harmonic_observables, anharmonic_observables/* or &
     &anharmonic_observables/*/temperature_* directory. The -d flag may be &
     &useful for this.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_dos_and_dispersion_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_dos_and_dispersion_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_dos_and_dispersion.py'),python_path)
end subroutine
end module
