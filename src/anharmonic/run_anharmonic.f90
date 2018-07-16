! ======================================================================
! Runs anharmonic DFT calculations, as set up by setup_anharmonic.
! ======================================================================
module run_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use anharmonic_common_module
  use polynomial_module
  
  use setup_anharmonic_module
  implicit none
  
  private
  
  public :: run_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_anharmonic() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'run_anharmonic'
  output%description = 'Runs DFT calculations set up by setup_anharmonic. &
     &Should be run after setup_anharmonic.'
  output%keywords = [                                                         &
     & KeywordData( 'calculations_to_run',                                    &
     &              'calculations_to_run specifies the first and last &
     &calculations to run inclusive. These should be given as two integers &
     &separated by a space. If not set then all calculations will be run. See &
     &calculation_directories.dat for the entire list of calculations to be &
     &run.',                                                                  &
     &              is_optional=.true.),                                      &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files.',                   &
     &               is_path=.true.),                                         &
     & KeywordData( 'no_cores',                                               &
     &              'no_cores is the number of cores on which DFT will be &
     &run. This is passed to the specified run script.',                      &
     &               default_value='1'),                                      &
     & KeywordData( 'calculation_type',                                       &
     &              'calculation_type specifies whether any electronic &
     &structure calculations should be run in addition to the user-defined &
     &script. Settings are: "none" and "quip".',                              &
     &              default_value='none') ]
  output%main_subroutine => run_anharmonic_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_anharmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory.
  type(String) :: wd
  
  ! Inputs
  integer, allocatable :: calculations_to_run(:)
  type(String)         :: run_script
  integer              :: no_cores
  type(String)         :: calculation_type
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: harmonic_path
  type(String)     :: file_type
  type(String)     :: seedname
  
  ! Calculation directories information.
  type(IFile)               :: calculation_directories_file
  type(String), allocatable :: calculation_directories(:)
  
  ! Electronic structure calculation runner.
  type(CalculationRunner) :: calculation_runner
  
  ! Temporary variables.
  integer :: i
  
  ! Read in inputs.
  wd = arguments%value('working_directory')
  if (arguments%is_set('calculations_to_run')) then
    calculations_to_run = &
       & int(split_line(arguments%value('calculations_to_run')))
  endif
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  calculation_type = arguments%value('calculation_type')
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic())
  call setup_anharmonic_arguments%read_file( &
     & wd//'/setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in calculation directories.
  calculation_directories_file = IFile(wd//'/calculation_directories.dat')
  calculation_directories = calculation_directories_file%lines()
  
  ! Select only those calculations specified by calculations_to_run,
  !    if relevant.
  if (allocated(calculations_to_run)) then
    if (size(calculations_to_run)/=2) then
      call print_line(ERROR//': Unable to parse calculations_to_run.')
      call print_line('calculations_to_run = "'//calculations_to_run//'"')
      stop
    elseif (calculations_to_run(1)<1) then
      call print_line(ERROR//': Index of first calculation to run less than &
         &0.')
      stop
    elseif (calculations_to_run(2)>size(calculation_directories)) then
      call print_line(ERROR//': Index of last calculation to run greater than &
         &no. calculations. See calculation_directories.dat for the entire &
         &list of calculations to be run.')
      stop
    elseif (calculations_to_run(2)<calculations_to_run(1)) then
      call print_line(ERROR//': Index of last calculation to run less than &
         &index of first calculation to run.')
      stop
    endif
    
    calculation_directories = calculation_directories(calculations_to_run(1): &
                                                     &calculations_to_run(2))
  endif
  
  ! Initialise calculation runner.
  calculation_runner = CalculationRunner(   &
     & working_directory = wd,              &
     & file_type         = file_type,       &
     & seedname          = seedname,        &
     & run_script        = run_script,      &
     & no_cores          = no_cores,        &
     & calculation_type  = calculation_type )
  
  ! Run calculations.
  call calculation_runner%run_calculations(calculation_directories)
end subroutine
end module
