! ======================================================================
! Plots the results of calculate_normal_modes.
! ======================================================================
module plot_normal_modes_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_normal_modes_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function plot_normal_modes_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_normal_modes'
  output%description = 'Plots the output of calculate_normal_modes. Should be &
     &run from within a qpoint_ directory. The -d flag may be useful for this.'
  output%keywords = plot_normal_modes_keywords()
  output%main_subroutine => plot_normal_modes
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_normal_modes(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  wd = arguments%value('working_directory')
  
  call execute_python(wd,str('plot_normal_modes.py'))
end subroutine
end module
