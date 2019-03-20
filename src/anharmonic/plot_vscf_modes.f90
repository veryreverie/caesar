! ======================================================================
! Plots the VSCF potential along each mode.
! ======================================================================
module plot_vscf_modes_module
  use common_module
  implicit none
  
  private
  
  public :: startup_plot_vscf_modes
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_vscf_modes()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_vscf_modes'
  mode%description = 'Plots the mapping of the vscf modes produced by &
     &map_vscf_modes.'
  mode%keywords = [                                                        &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_vscf_modes_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_vscf_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_vscf_modes.py'),python_path)
end subroutine
end module
