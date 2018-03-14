! ======================================================================
! Plots the results of calculate_normal_modes.
! ======================================================================
module plot_normal_modes_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_normal_modes_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                             &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
end function

function plot_normal_modes_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_normal_modes'
  output%description = 'Plots the output of calculate_normal_modes. Should be &
     &run from within a qpoint_ directory. The -d flag may be useful for this.'
  output%keywords = plot_normal_modes_keywords()
  output%main_subroutine => plot_normal_modes
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_normal_modes(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_normal_modes.py'),python_path)
end subroutine
end module
