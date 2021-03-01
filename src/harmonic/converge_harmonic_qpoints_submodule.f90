submodule (caesar_converge_harmonic_qpoints_module) caesar_converge_harmonic_qpoints_submodule
  use caesar_harmonic_module
contains

module procedure startup_converge_harmonic_qpoints
  type(CaesarMode) :: mode
  
  type(CaesarMode) :: setup_harmonic_mode
  type(CaesarMode) :: run_harmonic_mode
  type(CaesarMode) :: calculate_normal_modes_mode
  type(CaesarMode) :: calculate_harmonic_observables_mode
  
  setup_harmonic_mode = CaesarMode('setup_harmonic')
  run_harmonic_mode = CaesarMode('run_harmonic')
  calculate_normal_modes_mode = CaesarMode('calculate_normal_modes')
  calculate_harmonic_observables_mode = &
     & CaesarMode('calculate_harmonic_observables')
  
  mode%mode_name = 'converge_harmonic_qpoints'
  mode%description = 'Converges harmonic free energies w/r/t q-point spacing.'
  mode%keywords = [                                                           &
     & KeywordData( 'maximum_q-point_spacing',                                &
     &              'maximum_q-point_spacing is the largest q-point spacing &
     &at which energies will be calculated. maximum_q-point_spacing should be &
     &given in inverse Bohr.'),                                               &
     & KeywordData( 'no_q-point_spacings',                                    &
     &              'no_q-point_spacings is the approximate number of q-point &
     &spacings at which energies will be calculated. N.B. only q-point &
     &spacings corresponding to integer q-point grids will be sampled.'),     &
     & KeywordData( 'minimum_q-point_spacing',                                &
     &              'minimum_q-point_spacing is the smallest q-point spacing &
     &at which energies will be calculated. minimum_q-point_spacing should be &
     &given in inverse Bohr.'),                                               &
     & KeywordData( 'energy_tolerance',                                       &
     &              'energy_tolerance is the tolerance used in the &
     &convergence of the vibrational free energies. It is given in Hartree &
     &per primitive unit cell.',                                              &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'convergence_count',                                      &
     &              'convergence_count is the number of consecutive &
     &frequencies and/or energies that must be within the defined tolerance &
     &of each other for convergence to be considered reached',                &
     &              default_value='3')                                        ]
  mode%keywords = [ mode%keywords,                               &
                  & setup_harmonic_mode%keywords,                &
                  & run_harmonic_mode%keywords,                  &
                  & calculate_normal_modes_mode%keywords,        &
                  & calculate_harmonic_observables_mode%keywords ]
  mode%main_subroutine => converge_harmonic_qpoints_subroutine
  
  call mode%remove_keyword('q-point_grid')
  call mode%remove_keyword('supercells_to_run')
  call mode%remove_keyword('calculations_to_run')
  call mode%remove_keyword('exit_on_error')
  
  call add_mode(mode)
end procedure

module procedure converge_harmonic_qpoints_subroutine
  ! User inputs.
  type(String) :: seedname
  type(String) :: file_type
  real(dp)     :: maximum_spacing
  integer      :: no_spacings
  real(dp)     :: minimum_spacing
  logical      :: repeat_calculations
  real(dp)     :: energy_tolerance
  integer      :: convergence_count
  real(dp)     :: min_temperature
  real(dp)     :: max_temperature
  integer      :: random_seed
  
  ! Random generator for generating seed if not set.
  type(RandomReal) :: random_generator
  
  ! Structure data.
  type(String)        :: input_filename
  type(StructureData) :: structure
  
  ! q-point spacings and grids.
  real(dp),         allocatable :: qpoint_spacings(:)
  type(KpointGrid), allocatable :: qpoint_grids(:)
  logical,          allocatable :: unique_grids(:)
  
  ! Free energies.
  type(RealVector), allocatable :: free_energies(:)
  
  ! Directories and files.
  type(String) :: directory
  
  ! Temporary variables.
  logical :: converged
  integer :: i,ialloc
  
  ! Parse user inputs.
  seedname = arguments%value('seedname')
  file_type = arguments%value('file_type')
  maximum_spacing = dble(arguments%value('maximum_q-point_spacing'))
  no_spacings = int(arguments%value('no_q-point_spacings'))
  minimum_spacing = dble(arguments%value('minimum_q-point_spacing'))
  repeat_calculations = lgcl(arguments%value('repeat_calculations'))
  energy_tolerance = dble(arguments%value('energy_tolerance'))
  convergence_count = int(arguments%value('convergence_count'))
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  random_seed = random_generator%get_seed()
  
  ! Check user inputs.
  if (maximum_spacing<=minimum_spacing) then
    call print_line(ERROR//': maximum_spacing is smaller than &
       &minimum_spacing.')
    call quit()
  elseif (convergence_count<1) then
    call print_line(ERROR//': convergence_count must be at least 1.')
    call quit()
  endif
  
  ! Read in structure.
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Generate q-point grids.
  allocate( qpoint_spacings(no_spacings), &
          & qpoint_grids(no_spacings),    &
          & stat=ialloc); call err(ialloc)
  do i=1,no_spacings
    ! Evenly space q-point spacings in 1/spacing.
    qpoint_spacings(i) = (no_spacings-1)                   &
                     & / ( (no_spacings-i)/maximum_spacing &
                     &   + (i-1)/minimum_spacing           )
    
    ! Calculate the q-point grid corresponding to each spacing.
    qpoint_grids(i) = calculate_kpoint_grid( qpoint_spacings(i),           &
                                           & dble(structure%recip_lattice) )
    
    ! Re-calculate spacing to reflect the actual q-point grid.
    qpoint_spacings(i) = calculate_kpoint_spacing( &
                   & qpoint_grids(i),              &
                   & dble(structure%recip_lattice) )
  enddo
  
  ! Remove duplicate q-point grids.
  unique_grids = [ .true.,                                 &
                 & ( qpoint_grids(i)/=qpoint_grids(i-1),   &
                 &   i=2,                                  &
                 &   no_spacings                         ) ]
  qpoint_spacings = qpoint_spacings(filter(unique_grids))
  qpoint_grids = qpoint_grids(filter(unique_grids))
  no_spacings = size(qpoint_spacings)
  
  ! Run calculations.
  allocate(free_energies(0), stat=ialloc); call err(ialloc)
  converged = .false.
  do i=1,no_spacings
    directory = 'qpoints_'//join(str(qpoint_grids(i)%grid), delimiter='_')
    call mkdir(directory)
    
    free_energies = [ free_energies,                                &
                    & calculate_free_energies( directory,           &
                    &                          qpoint_grids(i),     &
                    &                          input_filename,      &
                    &                          repeat_calculations, &
                    &                          random_seed,         &
                    &                          arguments            ) ]
    
    if (i>convergence_count) then
      converged = all(                                                    &
         &   maximum_difference( free_energies(i),                        &
         &                       free_energies(i-convergence_count:i-1) ) &
         & < energy_tolerance                                             )
    endif
    
    call write_output_file( structure,        &
                          & energy_tolerance, &
                          & min_temperature,  &
                          & max_temperature,  &
                          & qpoint_spacings,  &
                          & free_energies     )
    
    call print_line('q-point spacing: '//qpoint_spacings(i)//' (Bohr^-1)')
    if (converged) then
      exit
    endif
  enddo
  
  if (converged) then
    call print_line('Convergence reached.')
  else
    call print_line('Convergence not reached.')
  endif
end procedure

module procedure calculate_free_energies
  type(String) :: output_file_name
  
  type(Dictionary) :: setup_harmonic_arguments
  type(Dictionary) :: run_harmonic_arguments
  type(Dictionary) :: calculate_normal_modes_arguments
  type(Dictionary) :: calculate_harmonic_observables_arguments
  
  type(IFile)                          :: thermodynamics_file
  type(String),            allocatable :: lines(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)
  
  ! Run Caesar.
  call copy_file(input_file_name, directory//'/'//input_file_name)
  output_file_name = &
     & directory//'/harmonic_observables/thermodynamic_variables.dat'
  if (repeat_calculations .or. .not. file_exists(output_file_name)) then
    setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
    call setup_harmonic_arguments%set(arguments)
    call setup_harmonic_arguments%set('working_directory', directory)
    call setup_harmonic_arguments%set( 'output_file', &
                                     & directory//'/setup_harmonic.output')
    call setup_harmonic_arguments%set('q-point_grid', str(qpoint_grid))
    call call_caesar('setup_harmonic', setup_harmonic_arguments)
    
    run_harmonic_arguments = Dictionary(CaesarMode('run_harmonic'))
    call run_harmonic_arguments%set(arguments)
    call run_harmonic_arguments%set('working_directory', directory)
    call run_harmonic_arguments%set( 'output_file', &
                                   & directory//'/run_harmonic.output')
    call call_caesar('run_harmonic', run_harmonic_arguments)
    
    calculate_normal_modes_arguments = &
       & Dictionary(CaesarMode('calculate_normal_modes'))
    call calculate_normal_modes_arguments%set(arguments)
    call calculate_normal_modes_arguments%set('working_directory', directory)
    call calculate_normal_modes_arguments%set(       &
       & 'output_file',                              &
       & directory//'/calculate_normal_modes.output' )
    call call_caesar( 'calculate_normal_modes',        &
                    & calculate_normal_modes_arguments )
    
    calculate_harmonic_observables_arguments = &
       & Dictionary(CaesarMode('calculate_harmonic_observables'))
    call calculate_harmonic_observables_arguments%set(arguments)
    call calculate_harmonic_observables_arguments%set( 'working_directory', &
                                                     & directory            )
    call calculate_harmonic_observables_arguments%set( 'random_seed',   &
                                                     & str(random_seed) )
    call calculate_harmonic_observables_arguments%set(       &
       & 'output_file',                                      &
       & directory//'/calculate_harmonic_observables.output' )
    call call_caesar( 'calculate_harmonic_observables',        &
                    & calculate_harmonic_observables_arguments )
  endif
  
  ! Read in thermodynamic data.
  thermodynamics_file = IFile(output_file_name)
  lines = thermodynamics_file%lines()
  thermodynamics = ThermodynamicData(lines(2:))
  
  output = vec(thermodynamics%free_energy)
end procedure

module procedure copy_file
  type(IFile) :: in_file
  type(OFile) :: out_file
  
  in_file = IFile(input)
  out_file = OFile(output)
  
  call out_file%print_lines(in_file%lines())
end procedure

module procedure write_output_file
  type(OFile) :: output_file
  
  integer :: i
  
  output_file = OFile('convergence.dat')
  
  call output_file%print_line('No. atoms                  : '// &
                             &structure%no_atoms)
  call output_file%print_line('Cell volume                : '//structure%volume)
  call output_file%print_line('Convergence threshold (Ha) : '//energy_tolerance)
  call output_file%print_line('Temperature range     (Ha) : '// &
                             &min_temperature//' '//max_temperature)
  call output_file%print_line('')
  
  call output_file%print_line('q-point spacing (Bohr^-1) | Free energies &
     &(Ha per primitive cell)')
  do i=1,size(free_energies)
    call output_file%print_line( qpoint_spacings(i)//' '// &
                               & free_energies(i)          )
  enddo
end procedure

module procedure maximum_difference
  if (size(this)/=size(that)) then
    call print_line(CODE_ERROR//': Vectors of different lengths.')
    call err()
  endif
  
  output = maxval(abs(dble(this-that)))
end procedure
end submodule
