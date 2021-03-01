submodule (caesar_run_anharmonic_module) caesar_run_anharmonic_submodule
  use caesar_anharmonic_module
contains

module procedure startup_run_anharmonic
  type(CaesarMode) :: mode
  
  mode%mode_name = 'run_anharmonic'
  mode%description = 'Runs DFT calculations set up by setup_anharmonic. &
     &Should be run after setup_anharmonic.'
  mode%keywords = [                                                           &
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
     &              'no_cores is the number of cores on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                       &
     & KeywordData( 'no_nodes',                                               &
     &              'no_nodes is the number of nodes on which the electronic &
     &structure calculation will be run. This is passed to the specified run &
     &script.',                                                               &
     &              default_value='1'),                                       &
     & KeywordData( 'run_script_data',                                        &
     &              'run_script_data will be passed to the specified run &
     &script after all other arguments. This should be used to pass &
     &information not covered by the other arguments.',                       &
     &              default_value=''),                                        &
     & KeywordData( 'calculation_type',                                       &
     &              'calculation_type specifies whether any electronic &
     &structure calculations should be run in addition to the user-defined &
     &script. Settings are: "none" and "quip".',                              &
     &              default_value='none'),                                    &
     & KeywordData( 'exit_on_error',                                          &
     &              'exit_on_error specifies whether or not the code will &
     &exit if the electronic structure run script returns a result code other &
     &than 0.',                                                               &
     &              default_value='false'),                                   &
     & KeywordData( 'repeat_calculations',                                    &
     &              'repeat_calculations specifies whether or not electronic &
     &calculations will be re-run if an electronic_structure.dat file is &
     &found in their directory.',                                             &
     &               default_value='true')                                    ]
  mode%main_subroutine => run_anharmonic_subroutine
  
  call add_mode(mode)
end procedure

module procedure run_anharmonic_subroutine
  ! Inputs
  integer, allocatable :: calculations_to_run(:)
  type(String)         :: run_script
  integer              :: no_cores
  integer              :: no_nodes
  type(String)         :: run_script_data
  type(String)         :: calculation_type
  logical              :: exit_on_error
  logical              :: repeat_calculations
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: file_type
  type(String)     :: seedname
  
  type(Dictionary) :: setup_anharmonic_arguments
  logical          :: use_forces
  logical          :: use_hessians
  logical          :: calculate_stress
  
  ! Calculation directories information.
  type(IFile)               :: calculation_directories_file
  type(String), allocatable :: calculation_directories(:)
  
  ! Electronic structure calculation runner.
  type(CalculationRunner) :: calculation_runner
  
  ! Read in inputs.
  if (arguments%is_set('calculations_to_run')) then
    calculations_to_run = &
       & int(split_line(arguments%value('calculations_to_run')))
  endif
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script_data = arguments%value('run_script_data')
  calculation_type = arguments%value('calculation_type')
  exit_on_error = lgcl(arguments%value('exit_on_error'))
  repeat_calculations = lgcl(arguments%value('repeat_calculations'))
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file( &
          & 'setup_harmonic.used_settings' )
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = Dictionary(CaesarMode('setup_anharmonic'))
  call setup_anharmonic_arguments%read_file('setup_anharmonic.used_settings')
  use_forces = lgcl(setup_anharmonic_arguments%value('use_forces'))
  use_hessians = lgcl(setup_anharmonic_arguments%value('use_hessians'))
  calculate_stress = lgcl(setup_anharmonic_arguments%value('calculate_stress'))
  
  ! Read in calculation directories.
  calculation_directories_file = IFile('calculation_directories.dat')
  calculation_directories = calculation_directories_file%lines()
  
  ! Select only those calculations specified by calculations_to_run,
  !    if relevant.
  if (allocated(calculations_to_run)) then
    if (size(calculations_to_run)/=2) then
      call print_line(ERROR//': Unable to parse calculations_to_run.')
      call print_line('calculations_to_run = "'//calculations_to_run//'"')
      call quit()
    elseif (calculations_to_run(1)<1) then
      call print_line(ERROR//': Index of first calculation to run less than &
         &0.')
      call quit()
    elseif (calculations_to_run(2)>size(calculation_directories)) then
      call print_line(ERROR//': Index of last calculation to run greater than &
         &no. calculations. See calculation_directories.dat for the entire &
         &list of calculations to be run.')
      call quit()
    elseif (calculations_to_run(2)<calculations_to_run(1)) then
      call print_line(ERROR//': Index of last calculation to run less than &
         &index of first calculation to run.')
      call quit()
    endif
    
    calculation_directories = calculation_directories( &
       & calculations_to_run(1):calculations_to_run(2) )
  endif
  
  ! Initialise calculation runner.
  calculation_runner = CalculationRunner(        &
     & file_type           = file_type,          &
     & seedname            = seedname,           &
     & run_script          = run_script,         &
     & no_cores            = no_cores,           &
     & no_nodes            = no_nodes,           &
     & run_script_data     = run_script_data,    &
     & calculation_type    = calculation_type,   &
     & use_forces          = use_forces,         &
     & use_hessians        = use_hessians,       &
     & calculate_stress    = calculate_stress,   &
     & exit_on_error       = exit_on_error,      &
     & repeat_calculations = repeat_calculations )
  
  ! Run calculations.
  call calculation_runner%run_calculations(calculation_directories)
end procedure
end submodule
