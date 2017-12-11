! ======================================================================
! Plots the results of converge_cutoff_and_kpoints.
! ======================================================================
module plot_cutoff_and_kpoints_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_cutoff_and_kpoints_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function plot_cutoff_and_kpoints_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_cutoff_and_kpoints'
  output%description = 'Plots the output of coverge_cutoff_and_kpoints'
  output%keywords = plot_cutoff_and_kpoints_keywords()
  output%main_subroutine => plot_cutoff_and_kpoints
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_cutoff_and_kpoints(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  wd = arguments%value('working_directory')
  
  call execute_python(wd,str('plot_cutoff_and_kpoints.py'))
end subroutine
end module