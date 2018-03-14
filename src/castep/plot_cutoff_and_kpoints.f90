! ======================================================================
! Plots the results of converge_cutoff_and_kpoints.
! ======================================================================
module plot_cutoff_and_kpoints_module
  use common_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_cutoff_and_kpoints_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                             &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
end function

function plot_cutoff_and_kpoints_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_cutoff_and_kpoints'
  output%description = 'Plots the output of coverge_cutoff_and_kpoints'
  output%keywords = plot_cutoff_and_kpoints_keywords()
  output%main_subroutine => plot_cutoff_and_kpoints
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_cutoff_and_kpoints(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_cutoff_and_kpoints.py'),python_path)
end subroutine
end module
