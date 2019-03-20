! ======================================================================
! Plots the anharmonic potential along multiple modes, along with the harmonic
!    and effective harmonic potentials along those modes mode..
! ======================================================================
module plot_potential_map_module
  use common_module
  implicit none
  
  private
  
  public :: startup_plot_potential_map
contains

! ----------------------------------------------------------------------
! Generates keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_plot_potential_map()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'plot_potential_map'
  mode%description = 'Plots the mapping of the anharmonic potential &
     &calculated by map_potential.'
  mode%keywords = [                                                      &
     & KeywordData( 'python_path',                                         &
     &              'python_path is the path to the Python 3 executable.', &
     &              default_value='python3') ]
  mode%main_subroutine => plot_potential_map_subroutine
  mode%suppress_settings_file = .true.
  
  call add_mode(mode)
end subroutine

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
