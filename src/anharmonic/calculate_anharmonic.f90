! ======================================================================
! Calculates anharmonic properties, using the results of run_anharmonic.
! ======================================================================
module calculate_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use shared_module
  use polynomial_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: calculate_anharmonic_keywords
  public :: calculate_anharmonic_mode
  public :: calculate_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function calculate_anharmonic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [ KeywordData:: ]
end function

function calculate_anharmonic_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'run_anharmonic'
  output%description = 'Uses the results of run_anharmonic to calculate &
     &anharmonic properties. Should be run after run_anharmonic.'
  output%keywords = calculate_anharmonic_keywords()
  output%main_subroutine => calculate_anharmonic
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: file_type
  type(String)     :: seedname
  real(dp)         :: symmetry_precision
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic_keywords())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
end subroutine
end module
