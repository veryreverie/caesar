! ======================================================================
! Plots the anharmonic potential along multiple modes, along with the harmonic
!    and effective harmonic potentials along those modes mode..
! ======================================================================
module plot_potential_map_module
  use common_module
  implicit none
  
  private
  
  public :: plot_potential_map
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
function plot_potential_map() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'plot_potential_map'
  output%description = 'Plots the mapping of the anharmonic potential &
     &calculated by map_potential.'
  output%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  output%main_subroutine => plot_potential_map_subroutine
  output%suppress_settings_file = .true.
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine plot_potential_map_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  type(String) :: python_path
  
  python_path = arguments%value('python_path')
  
  call execute_python(str('plot_potential_map.py'),python_path)
end subroutine
end module
