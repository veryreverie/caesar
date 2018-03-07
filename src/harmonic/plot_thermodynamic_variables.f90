! ======================================================================
! Plots the thermodynamic variables calculated by
!    calculate_harmonic_observables.
! ======================================================================
module plot_thermodynamic_variables_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_thermodynamic_variables_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [KeywordData::]
end function

function plot_thermodynamic_variables_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_thermodynamic_variables'
  output%description = 'Plots the thermodynamic variables &
     &calculated by calculate_harmonic_observables.'
  output%keywords = plot_thermodynamic_variables_keywords()
  output%main_subroutine => plot_thermodynamic_variables
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_thermodynamic_variables(arguments)
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: wd
  
  wd = arguments%value('working_directory')
  
  call execute_python(wd,str('plot_thermodynamic_variables.py'))
end subroutine
end module
