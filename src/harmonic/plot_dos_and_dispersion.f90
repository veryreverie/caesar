! ======================================================================
! Plots the results of calculate_dos_and_dispersion.
! ======================================================================
module plot_dos_and_dispersion_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_dos_and_dispersion_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function plot_dos_and_dispersion_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_dos_and_dispersion'
  output%description = 'Plots the output of calculate_dos_and_dispersion.'
  output%keywords = plot_dos_and_dispersion_keywords()
  output%main_subroutine => plot_dos_and_dispersion
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_dos_and_dispersion(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  wd = arguments%value('working_directory')
  
  call execute_python(wd,str('plot_dos_and_dispersion.py'))
end subroutine
end module
