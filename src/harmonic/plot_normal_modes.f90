! ======================================================================
! Plots the results of calculate_normal_modes.
! ======================================================================
module caesar_plot_normal_modes_module
  use caesar_common_module
  implicit none
  
  private
  
  public :: startup_plot_normal_modes
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_normal_modes()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_normal_modes'
  mode%description = 'Plots the output of calculate_normal_modes. Should be &
     &run from within a qpoint_ directory. The -d flag may be useful for this.'
  mode%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_normal_modes_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_normal_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_normal_modes.py'),python_path)
end subroutine
end module
