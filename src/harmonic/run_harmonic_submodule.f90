submodule (caesar_run_harmonic_module) caesar_run_harmonic_submodule
  use caesar_harmonic_module
contains

module procedure startup_run_harmonic
  type(CaesarMode) :: mode
  
  mode%mode_name = 'run_harmonic'
  mode%description = 'Runs DFT calculations set up by setup_harmonic. &
     &should be run after setup_harmonic.'
  mode%keywords = [                                                           &
     & KeywordData( 'supercells_to_run',                                      &
     &              'supercells_to_run is the indices of the first and last &
     &supercell to run. These should be specified as two integers separated &
     &by a space. If neither supercells_to_run nor directories_to_run is set, &
     &all calculations will be run.',                                         &
     &              exclusive_with=[str('calculations_to_run')]),             &
     & KeywordData( 'calculations_to_run', &
     &              'calculations_to_run is the indices of the first and last &
     &directories to run inclusive. These should be given as two integers &
     &separated by a space. If neither supercells_to_run nor &
     &directories_to_run is set, all calculations will be run. See &
     &harmonic_calculation_directories.dat for the entire list of &
     &calculations to be run.',                                               &
     &              exclusive_with=[str('supercells_to_run')]),               &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files.',                   &
     &              is_path=.true.),                                          &
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
  mode%main_subroutine => run_harmonic_subroutine
  
  call add_mode(mode)
end procedure

module procedure run_harmonic_subroutine
  ! User inputs.
  integer      :: supercells_to_run(2)
  integer      :: calculations_to_run(2)
  integer      :: no_cores
  integer      :: no_nodes
  type(String) :: run_script_data
  type(String) :: run_script
  type(String) :: calculation_type
  logical      :: exit_on_error
  logical      :: repeat_calculations
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  integer          :: no_supercells
  type(String)     :: file_type
  type(String)     :: seedname
  
  ! Electronic structure calculation runner.
  type(CalculationRunner) :: calculation_runner
  
  ! Atom and direction data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  integer                            :: atom
  type(String)                       :: direction
  type(String)                       :: atom_string
  
  ! Files and Directories.
  type(IFile)               :: no_supercells_file
  type(IFile)               :: unique_directions_file
  type(IFile)               :: calculation_directories_file
  type(String)              :: supercell_dir
  type(String)              :: run_dir
  type(String), allocatable :: run_dirs(:)
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script_data = arguments%value('run_script_data')
  calculation_type = arguments%value('calculation_type')
  exit_on_error = lgcl(arguments%value('exit_on_error'))
  repeat_calculations = lgcl(arguments%value('repeat_calculations'))
  
  ! --------------------------------------------------
  ! Read in arguments to previous calculations.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  no_supercells_file = IFile('no_supercells.dat')
  no_supercells = int(no_supercells_file%line(1))
  
  if (arguments%is_set('supercells_to_run')) then
    supercells_to_run = int(split_line(arguments%value('supercells_to_run')))
    calculations_to_run = [0,0]
  elseif (arguments%is_set('calculations_to_run')) then
    supercells_to_run = [0,0]
    calculations_to_run = int(split_line(       &
       & arguments%value('calculations_to_run') ))
  else
    supercells_to_run = [1, no_supercells]
    calculations_to_run = [0,0]
  endif
  
  ! --------------------------------------------------
  ! Check user inputs.
  ! --------------------------------------------------
  if (no_cores<=0) then
    call print_line('')
    call print_line('Error: no. cores must be >= 1')
    call err()
  elseif (no_nodes<=0) then
    call print_line('')
    call print_line('Error: no. nodes must be >= 1')
    call err()
  elseif (.not. file_exists(run_script)) then
    call print_line('')
    call print_line('Error: '//run_script//' does not exist.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation runner.
  ! --------------------------------------------------
  ! N.B. if Hessians are available, the alternative "read_normal_modes" should
  !    be used instead.
  calculation_runner = CalculationRunner(        &
     & file_type           = file_type,          &
     & seedname            = seedname,           &
     & run_script          = run_script,         &
     & no_cores            = no_cores,           &
     & no_nodes            = no_nodes,           &
     & run_script_data     = run_script_data,    &
     & calculation_type    = calculation_type,   &
     & use_forces          = .true.,             &
     & use_hessians        = .false.,            &
     & calculate_stress    = .false.,            &
     & exit_on_error       = exit_on_error,      &
     & repeat_calculations = repeat_calculations )
  
  ! --------------------------------------------------
  ! Calculate which directories to run.
  ! --------------------------------------------------
  if (.not. arguments%is_set('calculations_to_run')) then
    if (supercells_to_run(1)<=0) then
      call print_line('')
      call print_line('Error: first supercell must be > 0')
      call err()
    elseif (supercells_to_run(1)>no_supercells) then
      call print_line('')
      call print_line('Error: first supercell must be <= '//no_supercells)
      call err()
    elseif (supercells_to_run(2)<supercells_to_run(1)) then
      call print_line('')
      call print_line('Error: first supercell must be <= last supercell.')
      call err()
    elseif (supercells_to_run(2)>no_supercells) then
      call print_line('')
      call print_line('Error: last supercell must be <= '//no_supercells)
      call err()
    endif
    
    ! Loop over supercells.
    allocate(run_dirs(0), stat=ialloc); call err(ialloc)
    do i=supercells_to_run(1),supercells_to_run(2)
      supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
      
      unique_directions_file = IFile(supercell_dir//'/unique_directions.dat')
      unique_directions = UniqueDirection(unique_directions_file%sections())
      
      ! Loop over displacements within each supercell.
      do j=1,size(unique_directions)
        atom = unique_directions(j)%atom_id
        direction = unique_directions(j)%direction
        atom_string = left_pad(atom, str(maxval(unique_directions%atom_id)))
        
        run_dir = supercell_dir//'/atom.'//atom_string//'.'//direction
        
        ! Run calculation at each displacement.
        run_dirs = [run_dirs, run_dir]
      enddo
      
      deallocate(unique_directions, stat=ialloc); call err(ialloc)
    enddo
  else
    calculation_directories_file = IFile(       &
       & 'harmonic_calculation_directories.dat' )
    run_dirs = calculation_directories_file%lines()
    
    if (calculations_to_run(1)<0) then
      call print_line(ERROR//': first calculation to run must be >0.')
      call quit()
    elseif (calculations_to_run(2)>size(run_dirs)) then
      call print_line(ERROR//': last calculation to run must be <='// &
         & size(run_dirs)//'.')
      call quit()
    elseif (calculations_to_run(2)<calculations_to_run(1)) then
      call print_line(ERROR//': last calculation ro run must be >= first &
         &calculation to run.')
      call err()
    endif
    
    run_dirs = run_dirs(calculations_to_run(1):calculations_to_run(2))
  endif
  
  ! --------------------------------------------------
  ! Run calculations
  ! --------------------------------------------------
  call calculation_runner%run_calculations(run_dirs)
end procedure
end submodule
