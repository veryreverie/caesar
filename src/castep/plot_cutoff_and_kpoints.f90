! ======================================================================
! Plots the results of converge_cutoff_and_kpoints.
! ======================================================================
module plot_cutoff_and_kpoints_module
  use common_module
  implicit none
  
  private
  
  public :: plot_cutoff_and_kpoints
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_cutoff_and_kpoints() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_cutoff_and_kpoints'
  output%description = 'Plots the output of coverge_cutoff_and_kpoints'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_cutoff_and_kpoints_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_cutoff_and_kpoints_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  type(String) :: python_path
  
  wd = arguments%value('working_directory')
  python_path = arguments%value('python_path')
  
  call execute_python(wd,str('plot_cutoff_and_kpoints.py'),python_path)
end subroutine
end module
