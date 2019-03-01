! ======================================================================
! Converges harmonic frequencies and free energies w/r/t
!    CASTEP cutoff energy and k-point spacing.
! Runs multiple harmonic caesar calculations with different
!    cutoff energies and k-point spacings.
! ======================================================================
module converge_harmonic_frequencies_module
  use common_module

  implicit none
  
  private
  
  public :: converge_harmonic_frequencies
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext
! ----------------------------------------------------------------------
function converge_harmonic_frequencies() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'converge_harmonic_frequencies'
  output%description = 'Converges harmonic frequencies and free energies &
     &w/r/t cutoff energy, k-point spacing and electronic smearing. N.B. only &
     &the value being converged will be changed in each input file.'
  output%keywords = [                                                         &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Usual settings are: "castep", &
     &"caesar" and "xyz", but only CASTEP is supported by this convergence &
     &module at present',                                                     &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &              'seedname is the CASTEP seedname from which filenames &
     &are constructed.'),                                                     &
     & KeywordData( 'converge_cutoff',                                        &
     &              'converge_cutoff specifies whether or not cutoff energies &
     &will be converged.',                                                    &
     &              default_value='true'),                                    &
     & KeywordData( 'minimum_cutoff',                                         &
     &              'minimum_cutoff is the smallest cutoff energy which will &
     &be tested. minimum_cutoff should be given in Hartree.'),                &
     & KeywordData( 'no_cutoffs',                                             &
     &              'no_cutoffs is the number of cutoff energies which will &
     &be sampled.'),                                                          &
     & KeywordData( 'maximum_cutoff',                                         &
     &              'maximum_cutoff is the cutoff energy at which calculation &
     &will be terminated if convergence has not been reached. maximum_cutoff &
     &should be given in Hartree.'),                                          &
     & KeywordData( 'converge_k-point_spacing',                               &
     &              'converge_k-point_spacing specifies whether or not &
     &k-point spacing will be converged.',                                    &
     &              default_value='true'),                                    &
     & KeywordData( 'maximum_k-point_spacing',                                &
     &              'maximum_k-point_spacing is the largest k-point spacing &
     &which will be sampled. maximum_k-point_spacing should be given in &
     &inverse Bohr.'),                                                        &
     & KeywordData( 'no_k-point_spacings',                                    &
     &              'no_k-point_spacings is the approximate number of k-point &
     &spacings which will be sampled. N.B. since only k-point spacings &
     &corresponding to integer k-point grids will be sampled.'),              &
     & KeywordData( 'minimum_k-point_spacing',                                &
     &              'minimum_k-point_spacing is the k-point spacing at which &
     &calculation will be terminated if convergence has not been reached. &
     &minimum_k-point_spacing should be given in inverse Bohr.'),             &
     & KeywordData( 'k-point_parity',                                         &
     &              'k-point_parity specifies which k-point grids are &
     &sampled. If k-point_parity is "odd" or "even" then only k-point grids &
     &with an odd or even number of k-points in all directions are sampled. &
     &If k-point_parity is "both" then all k-point grids are samples.',       &
     &              default_value='both'),                                    &
     & KeywordData( 'converge_smearing',                                      &
     &              'converge_smearing specifies whether or not electron &
     &smearing will be converged.',                                           &
     &              default_value='true'),                                    &
     & KeywordData( 'maximum_smearing',                                       &
     &              'maximum_smearing is the largest electronic smearing &
     &which will be sampled. maximum_smearing should be given in Hartree.'),  &
     & KeywordData( 'no_smearings',                                           &
     &              'no_smearings is the number of electronic smearings which &
     &which will be sampled.'),                                               &
     & KeywordData( 'minimum_smearing',                                       &
     &              'minimum_smearing is the electronic smearing at which &
     &calculation will be terminated if convergence has not been reached. &
     &minimum_smearing should be given in Hartree.'),                         &
     & KeywordData( 'converge_energies',                                      &
     &              'converge_energies determines whether convergence of the &
     &free energies is monitored as well as convergence of harmonic &
     &frequencies.',                                                          &
     &              default_value='false'),                                   &
     & KeywordData( 'freq_tolerance',                                         &
     &              'freq_tolerance is the tolerance used in the convergence &
     &of the harmonic frequencies.',                                          &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'energy_tolerance',                                       &
     &              'energy_tolerance is the tolerance used in the &
     &convergence of the vibrational free energies. It is given in Hartree &
     &per primitive unit cell.',                                              &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'convergence_count',                                      &
     &              'convergence_count is the number of consecutive &
     &frequencies and/or energies that must be within the defined tolerance &
     &of each other for convergence to be considered reached',                &
     &              default_value='3'),                                       &
     & KeywordData( 'q-point_grid',                                           &
     &              'q-point_grid is the number of q-points in each direction &
     &in a Monkhorst-Pack grid. This should be specified as three integers &
     &separated by spaces.'),                                                 &
     & KeywordData( 'symmetry_precision',                                     &
     &              'In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. symmetry_precision should be given in Bohr.', &
     &              default_value='0.1'),                                     &
     & KeywordData( 'harmonic_displacement',                                  &
     &              'harmonic_displacement is the distance in bohr by which &
     &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
     &surface.',                                                              &
     &              default_value='0.01'),                                    &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files. This script should &
     &have executable file permissions and should run CASTEP',                &
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
     & KeywordData( 'repeat_calculations',                                    &
     &              'repeat_calculations specifies whether or not electronic &
     &calculations will be re-run if an electronic_structure.dat file is &
     &found in their directory.',                                             &
     &               default_value='true'),                                   &
     & KeywordData( 'acoustic_sum_rule',                                      &
     &              'acoustic_sum_rule specifies where the acoustic sum rule &
     &is applied. The options are "off", "forces", "matrices" and "both".',   &
     &              default_value='both'),                                    &
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                      &
     &              default_value='0'),                                       &
     & KeywordData( 'max_temperature',                                        &
     &              'max_temperature is the maximum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.',                                                      &
     &              default_value='500'),                                     &
     & KeywordData( 'no_temperature_steps',                                   &
     &              'no_temperature_steps is the number of temperatures at &
     &which thermodynamic quantities are calculated.',                        &
     &              default_value='6'),                                       &
     & KeywordData( 'min_frequency',                                          &
     &              'min_frequency is the frequency below which modes will be &
     &ignored when calculating thermodynamic quantities. min_frequency should &
     &be given in Hartree.',                                                  &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'path',                                                   &
     &              'path is the path through fractional reciprocal space &
     &which will be mapped by the phonon dispersion curve. The path should be &
     &specified as labels and q-points, separated by commas. The Gamma-point &
     &should be labelled G.',                                                 &
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5,              &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5'),                          &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='100000')                                   ]
  output%main_subroutine => converge_harmonic_frequencies_subroutine
end function

! ----------------------------------------------------------------------
! Main program
! ----------------------------------------------------------------------
subroutine converge_harmonic_frequencies_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs
  type(String) :: wd
  type(String) :: seedname
  type(String) :: file_type
  
  logical      :: converge_cutoff
  real(dp)     :: minimum_cutoff
  integer      :: no_cutoffs
  real(dp)     :: maximum_cutoff
  
  logical      :: converge_spacing
  real(dp)     :: maximum_spacing
  integer      :: no_spacings
  real(dp)     :: minimum_spacing
  type(String) :: kpoint_parity
  
  logical      :: converge_smearing
  real(dp)     :: minimum_smearing
  integer      :: no_smearings
  real(dp)     :: maximum_smearing
  
  integer      :: grid(3)
  real(dp)     :: symmetry_precision
  real(dp)     :: harmonic_displacement
  integer      :: no_cores
  integer      :: no_nodes
  type(String) :: run_script
  logical      :: converge_energies
  logical      :: repeat_calculations
  type(String) :: acoustic_sum_rule
  real(dp)     :: min_temperature
  real(dp)     :: max_temperature
  integer      :: no_temperature_steps
  real(dp)     :: min_frequency
  type(String) :: path
  integer      :: no_dos_samples
  real(dp)     :: freq_tolerance
  real(dp)     :: energy_tolerance
  integer      :: convergence_count
  integer      :: random_seed
  
  ! Structure.
  type(String)        :: input_filename
  type(StructureData) :: structure

  ! Energy cut-offs and kpoint spacings
  real(dp),         allocatable :: cutoffs(:)
  real(dp),         allocatable :: kpoint_spacings(:)
  type(KpointGrid), allocatable :: kpoint_grids(:)
  logical,          allocatable :: unique_grids(:)
  real(dp),         allocatable :: smearings(:)
  
  ! Variables to converge.
  type(RealVector), allocatable :: cutoff_frequencies(:)
  type(RealVector), allocatable :: cutoff_free_energies(:)
  type(RealVector), allocatable :: kpoint_frequencies(:)
  type(RealVector), allocatable :: kpoint_free_energies(:)
  type(RealVector), allocatable :: smearing_frequencies(:)
  type(RealVector), allocatable :: smearing_free_energies(:)
  
  ! Working variables
  integer :: no_qpoints
  logical :: frequencies_converged
  logical :: energies_converged
  
  ! Random generator for generating seed if not set.
  type(RandomReal) :: random_generator
  
  ! Files and directories.
  type(String) :: directory
  
  ! Temporary variables
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  
  converge_cutoff = lgcl(arguments%value('converge_cutoff'))
  minimum_cutoff = dble(arguments%value('minimum_cutoff'))
  no_cutoffs = int(arguments%value('no_cutoffs'))
  maximum_cutoff = dble(arguments%value('maximum_cutoff'))
  
  converge_spacing = lgcl(arguments%value('converge_k-point_spacing'))
  maximum_spacing = dble(arguments%value('maximum_k-point_spacing'))
  no_spacings = int(arguments%value('no_k-point_spacings'))
  minimum_spacing = dble(arguments%value('minimum_k-point_spacing'))
  kpoint_parity = arguments%value('k-point_parity')
  
  converge_smearing = lgcl(arguments%value('converge_smearing'))
  maximum_smearing = dble(arguments%value('maximum_smearing'))
  no_smearings = int(arguments%value('no_smearings'))
  minimum_smearing = dble(arguments%value('minimum_smearing'))
  
  grid = int(split_line(arguments%value('q-point_grid')))
  no_qpoints = grid(1)*grid(2)*grid(3)
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script = arguments%value('run_script')
  converge_energies = lgcl(arguments%value('converge_energies'))
  repeat_calculations = lgcl(arguments%value('repeat_calculations'))
  acoustic_sum_rule = arguments%value('acoustic_sum_rule')
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path = arguments%value('path')
  no_dos_samples = int(arguments%value('no_dos_samples'))
  freq_tolerance = dble(arguments%value('freq_tolerance'))
  energy_tolerance = dble(arguments%value('energy_tolerance'))
  convergence_count = int(arguments%value('convergence_count'))
  
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  random_seed = random_generator%get_seed()
  
  if (maximum_cutoff<=minimum_cutoff) then
    call print_line(ERROR//': maximum_cutoff is smaller than minimum_cutoff.')
    call quit()
  elseif (maximum_spacing<=minimum_spacing) then
    call print_line(ERROR//': maximum_k-point_spacing is smaller than &
       &minimum_k-point_spacing.')
    call quit()
  elseif (maximum_smearing<=minimum_smearing) then
    call print_line(ERROR//': maximum_smearing is smaller than &
       &minimum_smearing.')
    call quit()
  elseif (file_type/='castep' .and. file_type/='quantum_espresso') then
    call print_line(ERROR//': castep and quantum_espresso are the only &
       &accepted file types for this mode.')
    call quit()
  elseif (max_temperature<=min_temperature) then
    call print_line(ERROR//': max_temperature is smaller than &
       &min_temperature.')
    call quit()
  elseif (convergence_count<1) then
    call print_line(ERROR//': convergence_count must be at least 1.')
    call quit()
  endif
  
  ! Read in structure.
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! Initialise arrays.
  cutoff_frequencies = [RealVector::]
  cutoff_free_energies = [RealVector::]
  kpoint_frequencies = [RealVector::]
  kpoint_free_energies = [RealVector::]
  smearing_frequencies = [RealVector::]
  smearing_free_energies = [RealVector::]
  
  ! Calculate spacings, smearings and cutoffs.
  allocate( cutoffs(no_cutoffs),          &
          & kpoint_spacings(no_spacings), &
          & kpoint_grids(no_spacings),    &
          & smearings(no_smearings),      &
          & stat=ialloc); call err(ialloc)
  do i=1,no_cutoffs
    cutoffs(i) = ( (no_cutoffs-i)*minimum_cutoff   &
             &   + (i-1)*maximum_cutoff          ) &
             & / (no_cutoffs-1)
  enddo
  
  do i=1,no_smearings
    smearings(i) = ( (no_smearings-i)*maximum_smearing &
               &   + (i-1)*minimum_smearing          ) &
               & / (no_smearings-1)
  enddo
  
  do i=1,no_spacings
    ! Evenly space k-point spacings in 1/spacing co-ordinates.
    kpoint_spacings(i) = 1                                     &
                     & / ( ( (no_spacings-i)/maximum_spacing   &
                     &     + (i-1)/minimum_spacing           ) &
                     &   / (no_spacings-1)                     )
    
    ! Calculate the k-point grid corresponding to each spacing.
    kpoint_grids(i) = calculate_kpoint_grid( kpoint_spacings(i),           &
                                           & dble(structure%recip_lattice) )
    
    ! Increase any grid elements which have the wrong parity.
    if (kpoint_parity=='both') then
      continue
    elseif (kpoint_parity=='odd') then
      do j=1,3
        if (modulo(kpoint_grids(i)%grid(j),2)==0) then
          kpoint_grids(i)%grid(j) = kpoint_grids(i)%grid(j)+1
        endif
      enddo
    elseif (kpoint_parity=='even') then
      do j=1,3
        if (modulo(kpoint_grids(i)%grid(j),2)==1) then
          kpoint_grids(i)%grid(j) = kpoint_grids(i)%grid(j)+1
        endif
      enddo
    else
      call print_line(ERROR//': Unexpected value for k-point_parity.')
      call quit()
    endif
    
    ! Re-calculate spacing to reflect the actual k-point grid.
    kpoint_spacings(i) = calculate_kpoint_spacing( &
                   & kpoint_grids(i),              &
                   & dble(structure%recip_lattice) )
  enddo
  
  ! Remove duplicate k-point grids.
  unique_grids = [(.true., i=1, no_spacings)]
  do i=2,no_spacings
    unique_grids(i) = .not. all(kpoint_grids(i)%grid==kpoint_grids(i-1)%grid)
  enddo
  kpoint_spacings = kpoint_spacings(filter(unique_grids))
  kpoint_grids = kpoint_grids(filter(unique_grids))
  no_spacings = size(kpoint_spacings)
  
  ! --------------------------------------------------
  ! Run cut-off energy convergence.
  ! --------------------------------------------------
  if (converge_cutoff) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over cutoffs, running calculations at each.
    do i=1,no_cutoffs
      directory = 'cutoff_'//trim(str(                           &
         & cutoffs(i),                                           &
         & settings=PrintSettings( decimal_places        = 2,    &
         &                         floating_point_format = 'f' ) ))
      call mkdir(directory)
      
      cutoff_frequencies = [                                 &
         & cutoff_frequencies,                               &
         & calculate_frequencies(                            &
         &    directory             = directory,             &
         &    cutoff                = cutoffs(i),            &
         &    seedname              = seedname,              &
         &    no_qpoints            = no_qpoints,            &
         &    file_type             = file_type,             &
         &    grid                  = grid,                  &
         &    symmetry_precision    = symmetry_precision,    &
         &    harmonic_displacement = harmonic_displacement, &
         &    run_script            = run_script,            &
         &    no_cores              = no_cores,              &
         &    no_nodes              = no_nodes,              &
         &    repeat_calculations   = repeat_calculations,   &
         &    acoustic_sum_rule     = acoustic_sum_rule      ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         cutoff_frequencies(i),                        &
           &         cutoff_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        cutoff_free_energies = [                            &
           & cutoff_free_energies,                          &
           & calculate_free_energies( directory,            &
           &                          min_temperature,      &
           &                          max_temperature,      &
           &                          no_temperature_steps, &
           &                          min_frequency,        &
           &                          path,                 &
           &                          no_dos_samples,       &
           &                          random_seed,          &
           &                          repeat_calculations   ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         cutoff_free_energies(i),                        &
             &         cutoff_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      ! Check for convergence.
      call print_line('Cut-off energy: '//cutoffs(i)//' (Ha)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
  endif
  
  ! --------------------------------------------------
  ! Run k-points spacing convergence.
  ! --------------------------------------------------
  if (converge_spacing) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over spacings, running calculations at each.
    do i=1,no_spacings
      directory = 'kpoints_'//join(str(kpoint_grids(i)%grid), delimiter='_')
      call mkdir(directory)
      
      kpoint_frequencies = [                                 &
         & kpoint_frequencies,                               &
         & calculate_frequencies(                            &
         &    directory             = directory,             &
         &    kpoint_grid           = kpoint_grids(i),       &
         &    seedname              = seedname,              &
         &    no_qpoints            = no_qpoints,            &
         &    file_type             = file_type,             &
         &    grid                  = grid,                  &
         &    symmetry_precision    = symmetry_precision,    &
         &    harmonic_displacement = harmonic_displacement, &
         &    run_script            = run_script,            &
         &    no_cores              = no_cores,              &
         &    no_nodes              = no_nodes,              &
         &    repeat_calculations   = repeat_calculations,   &
         &    acoustic_sum_rule     = acoustic_sum_rule      ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         kpoint_frequencies(i),                        &
           &         kpoint_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        kpoint_free_energies = [                            &
           & kpoint_free_energies,                          &
           & calculate_free_energies( directory,            &
           &                          min_temperature,      &
           &                          max_temperature,      &
           &                          no_temperature_steps, &
           &                          min_frequency,        &
           &                          path,                 &
           &                          no_dos_samples,       &
           &                          random_seed,          &
           &                          repeat_calculations   ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         kpoint_free_energies(i),                        &
             &         kpoint_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      call print_line('K-point spacing: '//kpoint_spacings(i)//' (Bohr^-1)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
  endif
  
  ! Run electronic smearing convergence.
  if (converge_smearing) then
    ! Initialise convergence checks and observable arrays.
    frequencies_converged = .false.
    energies_converged = .false.
    if (.not. converge_energies) then
     energies_converged = .true.
    endif
    
    ! Loop over cutoffs, running calculations at each.
    do i=1,no_smearings
      directory = 'smearing_'//trim(str(                         &
         & smearings(i),                                         &
         & settings=PrintSettings( decimal_places        = 6,    &
         &                         floating_point_format = 'f' ) ))
      call mkdir(directory)
      
      smearing_frequencies = [                               &
         & smearing_frequencies,                             &
         & calculate_frequencies(                            &
         &    directory             = directory,             &
         &    smearing              = smearings(i),          &
         &    seedname              = seedname,              &
         &    no_qpoints            = no_qpoints,            &
         &    file_type             = file_type,             &
         &    grid                  = grid,                  &
         &    symmetry_precision    = symmetry_precision,    &
         &    harmonic_displacement = harmonic_displacement, &
         &    run_script            = run_script,            &
         &    no_cores              = no_cores,              &
         &    no_nodes              = no_nodes,              &
         &    repeat_calculations   = repeat_calculations,   &
         &    acoustic_sum_rule     = acoustic_sum_rule      ) ]
      
      if (i>convergence_count) then
        frequencies_converged =                                    &
           & all( maximum_difference(                              &
           &         smearing_frequencies(i),                        &
           &         smearing_frequencies(i-convergence_count:i-1) ) &
           &    < freq_tolerance                                   )
      endif

      if (converge_energies .and. .not. energies_converged) then
        smearing_free_energies = [                          &
           & smearing_free_energies,                        &
           & calculate_free_energies( directory,            &
           &                          min_temperature,      &
           &                          max_temperature,      &
           &                          no_temperature_steps, &
           &                          min_frequency,        &
           &                          path,                 &
           &                          no_dos_samples,       &
           &                          random_seed,          &
           &                          repeat_calculations   ) ]

        if (i>convergence_count) then
          energies_converged =                                         &
             & all( maximum_difference(                                &
             &         smearing_free_energies(i),                        &
             &         smearing_free_energies(i-convergence_count:i-1) ) &
             &    < energy_tolerance                                   )
        endif
      endif
      
      ! Write output file.
      call write_output_file( cutoffs,               &
                            & cutoff_frequencies,    &
                            & cutoff_free_energies,  &
                            & kpoint_spacings,       &
                            & kpoint_frequencies,    &
                            & kpoint_free_energies,  &
                            & smearings,             &
                            & smearing_frequencies,  &
                            & smearing_free_energies )
      
      ! Check for convergence.
      call print_line('Electronic smearing: '//smearings(i)//' (Ha)')
      if (frequencies_converged .and. energies_converged) then
        call print_line('Convergence reached.')
        exit
      endif
    enddo
    
    if (.not. (frequencies_converged .or. energies_converged)) then
      call print_line('Convergence not reached.')
    endif
    
  endif
end subroutine

subroutine write_output_file(cutoffs,cutoff_frequencies,cutoff_free_energies, &
   & kpoint_spacings,kpoint_frequencies,kpoint_free_energies,smearings,       &
   & smearing_frequencies,smearing_free_energies                              )
  implicit none
  
  real(dp),         intent(in) :: cutoffs(:)
  type(RealVector), intent(in) :: cutoff_frequencies(:)
  type(RealVector), intent(in) :: cutoff_free_energies(:)
  real(dp),         intent(in) :: kpoint_spacings(:)
  type(RealVector), intent(in) :: kpoint_frequencies(:)
  type(RealVector), intent(in) :: kpoint_free_energies(:)
  real(dp),         intent(in) :: smearings(:)
  type(RealVector), intent(in) :: smearing_frequencies(:)
  type(RealVector), intent(in) :: smearing_free_energies(:)
  
  type(OFile)  :: output_file
  
  integer :: i
  
  ! Open output file.
  output_file = OFile('convergence.dat')
  
  if (size(cutoff_frequencies)>0) then
    call output_file%print_line('Cutoff energy (Ha) | Mode frequencies (Ha)')
    do i=1,size(cutoff_frequencies)
      call output_file%print_line(cutoffs(i)//' '//cutoff_frequencies(i))
    enddo
  endif
  
  if (size(cutoff_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('Cutoff energy (Ha) | Free energies &
       &(Ha per primitive cell)')
    do i=1,size(cutoff_free_energies)
      call output_file%print_line(cutoffs(i)//' '//cutoff_free_energies(i))
    enddo
  endif
  
  if (size(cutoff_frequencies)>0 .and. size(kpoint_frequencies)>0) then
    call output_file%print_line('')
  endif
  
  if (size(kpoint_frequencies)>0) then
    call output_file%print_line('k-point spacing (Bohr^-1) | Mode frequencies &
       &(Ha)')
    do i=1,size(kpoint_frequencies)
      call output_file%print_line( kpoint_spacings(i) //' '// &
                                 & kpoint_frequencies(i)      )
    enddo
  endif
  
  if (size(kpoint_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('k-point spacing (Bohr^-1) | Free energies &
       &(Ha per primitive cell)')
    do i=1,size(kpoint_free_energies)
      call output_file%print_line( kpoint_spacings(i)//' '// &
                                 & kpoint_free_energies(i)   )
    enddo
  endif
  
  if ( ( size(cutoff_frequencies)>0 .or.         &
     &   size(kpoint_frequencies)>0      ) .and. &
     & size(smearing_frequencies)>0              ) then
    call output_file%print_line('')
  endif
  
  if (size(smearing_frequencies)>0) then
    call output_file%print_line('Electronic smearing (Hartree) | Mode &
       &frequencies (Ha)')
    do i=1,size(smearing_frequencies)
      call output_file%print_line( smearings(i)    //' '// &
                                 & smearing_frequencies(i) )
    enddo
  endif
  
  if (size(smearing_free_energies)>0) then
    call output_file%print_line('')
    call output_file%print_line('Electronic smearing (Hartree) | Free &
       &energies (Ha per primitive cell)')
    do i=1,size(smearing_free_energies)
      call output_file%print_line( smearings(i)      //' '// &
                                 & smearing_free_energies(i) )
    enddo
  endif
end subroutine

function calculate_frequencies(directory,cutoff,kpoint_grid,smearing,        &
   & seedname,no_qpoints,file_type,grid,symmetry_precision,                  &
   & harmonic_displacement,run_script,no_cores,no_nodes,repeat_calculations, &
   & acoustic_sum_rule) result(output)
  implicit none
  
  type(String),        intent(in)           :: directory
  real(dp),            intent(in), optional :: cutoff
  type(KpointGrid),    intent(in), optional :: kpoint_grid
  real(dp),            intent(in), optional :: smearing
  type(String),        intent(in)           :: seedname
  integer,             intent(in)           :: no_qpoints
  type(String),        intent(in)           :: file_type
  integer,             intent(in)           :: grid(3)
  real(dp),            intent(in)           :: symmetry_precision
  real(dp),            intent(in)           :: harmonic_displacement
  type(String),        intent(in)           :: run_script
  integer,             intent(in)           :: no_cores
  integer,             intent(in)           :: no_nodes
  logical,             intent(in)           :: repeat_calculations
  type(String),        intent(in)           :: acoustic_sum_rule
  type(RealVector)                          :: output
  
  type(String) :: qpoint_dir
  type(IFile)  :: complex_modes_file
  
  type(ComplexMode), allocatable :: modes(:)
  
  integer :: i
  
  ! Check input.
  if (count([ present(cutoff),      &
            & present(kpoint_grid), &
            & present(smearing)     ])/=1) then
    call print_line(CODE_ERROR//': calculate_frequencies must be passed &
       &exactly one of cutoff, kpoint_grid and smearing.')
  endif
  
  ! Write DFT files
  if (file_type=='castep') then
    call write_castep_files(directory, seedname, cutoff, kpoint_grid, smearing)
  elseif (file_type=='quantum_espresso') then
    call write_qe_file(directory, seedname, cutoff, kpoint_grid, smearing)
  endif
  
  ! Call setup_harmonic.
  call call_caesar(                                              &
     & 'setup_harmonic -d '//directory                   //' '// &
     & '--file_type '//file_type                         //' '// &
     & '--seedname '//seedname                           //' '// &
     & '--q-point_grid '//grid                           //' '// &
     & '--symmetry_precision '//symmetry_precision       //' '// &
     & '--harmonic_displacement '//harmonic_displacement         )
  
  ! Call run_harmonic.
  call call_caesar( 'run_harmonic -d '//directory         //' '// &
                  & '--run_script '//run_script           //' '// &
                  & '--no_cores '//no_cores               //' '// &
                  & '--no_nodes '//no_nodes               //' '// &
                  & '--exit_on_error true'                //' '// &
                  & '--repeat_calculations '//repeat_calculations )
  
  ! Call calculate_normal_modes.
  call call_caesar( 'calculate_normal_modes -d '//directory   //' '// &
                  & '--acoustic_sum_rule '//acoustic_sum_rule         )
  
  ! Read in normal modes.
  modes = [ComplexMode::]
  do i=1,no_qpoints
    qpoint_dir = directory//'/qpoint_'//left_pad(i, str(no_qpoints))
    complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    modes = [modes, ComplexMode(complex_modes_file%sections())]
  enddo
  
  output = vec(modes%frequency)
end function

function calculate_free_energies(directory,min_temperature,max_temperature, &
   & no_temperature_steps,min_frequency,path,no_dos_samples,random_seed,    &
   & repeat_calculations) result(output)
  implicit none
  
  type(String), intent(in) :: directory
  real(dp),     intent(in) :: min_temperature
  real(dp),     intent(in) :: max_temperature
  integer,      intent(in) :: no_temperature_steps
  real(dp),     intent(in) :: min_frequency
  type(String), intent(in) :: path
  integer,      intent(in) :: no_dos_samples
  integer,      intent(in) :: random_seed
  logical,      intent(in) :: repeat_calculations
  type(RealVector)         :: output
  
  type(String)                         :: file_name
  type(IFile)                          :: thermodynamics_file
  type(String),            allocatable :: lines(:)
  type(ThermodynamicData), allocatable :: thermodynamics(:)

  ! Call calculate_harmonic_observables.
  file_name = directory//'/harmonic_observables/thermodynamic_variables.dat'
  if (repeat_calculations .or. .not. file_exists(file_name)) then
    call call_caesar(                                            &
       & 'calculate_harmonic_observables -d '//directory //' '// &
       & '--min_temperature '//min_temperature           //' '// &
       & '--max_temperature '//max_temperature           //' '// &
       & '--no_temperature_steps '//no_temperature_steps //' '// &
       & '--min_frequency '//min_frequency               //' '// &
       & '--path '//path                                 //' '// &
       & '--no_dos_samples '//no_dos_samples             //' '// &
       & '--random_seed '//random_seed                           )
  endif
  
  ! Read in thermodynamic data.
  thermodynamics_file = IFile(file_name)
  lines = thermodynamics_file%lines()
  thermodynamics = ThermodynamicData(lines(2:))
  
  output = vec(thermodynamics%free_energy)
end function

subroutine write_castep_files(directory,seedname,cutoff,kpoint_grid,smearing)
  implicit none
  
  type(String),     intent(in)           :: directory
  type(String),     intent(in)           :: seedname
  real(dp),         intent(in), optional :: cutoff
  type(KpointGrid), intent(in), optional :: kpoint_grid
  real(dp),         intent(in), optional :: smearing
  
  type(IFile) :: input_cell_file
  type(OFile) :: output_cell_file
  type(IFile) :: input_param_file
  type(OFile) :: output_param_file
  
  logical :: line_replaced
  
  type(String), allocatable :: lines(:)
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! --------------------------------------------------
  ! Update cell file.
  ! --------------------------------------------------
  
  ! Read input cell file.
  input_cell_file = IFile(seedname//'.cell')
  lines = input_cell_file%lines()
  
  ! Update k-point spacing.
  if (present(kpoint_grid)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (size(line)>=1) then
        if ( line(1)=='kpoint_mp_spacing'  .or. &
           & line(1)=='kpoint_mp_grid'     .or. &
           & line(1)=='kpoints_mp_spacing' .or. &
           & line(1)=='kpoints_mp_grid'         ) then
          lines(i) = 'kpoint_mp_grid '//kpoint_grid
          line_replaced = .true.
          exit
        endif
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'kpoint_mp_grid '//kpoint_grid]
    endif
  endif
  
  ! Write output cell file.
  output_cell_file = OFile(directory//'/'//seedname//'.cell')
  call output_cell_file%print_lines(lines)
  
  ! --------------------------------------------------
  ! Update param file.
  ! --------------------------------------------------
  
  ! Read input param file.
  input_param_file = IFile(seedname//'.param')
  lines = input_param_file%lines()
  
  ! Update electronic cut-off energy.
  if (present(cutoff)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (line(1)=='cut_off_energy') then
        lines(i) = 'cut_off_energy '//cutoff*EV_PER_HARTREE
        line_replaced = .true.
        exit
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'cut_off_energy '//cutoff*EV_PER_HARTREE]
    endif
  endif
  
  ! Update electronic smearing.
  if (present(smearing)) then
    line_replaced = .false.
    do i=1,size(lines)
      line = split_line(lower_case(lines(i)))
      if (line(1)=='smearing_width') then
        lines(i) = 'smearing_width '//smearing/KB_IN_AU//' K'
        line_replaced = .true.
        exit
      endif
    enddo
    
    if (.not. line_replaced) then
      lines = [lines, 'smearing_width '//smearing/KB_IN_AU//' K']
    endif
  endif
  
  ! Write output param file.
  output_param_file = OFile(directory//'/'//seedname//'.param')
  call output_param_file%print_lines(lines)
end subroutine

subroutine write_qe_file(directory,seedname,cutoff,kpoint_grid,smearing)
  implicit none
  
  type(String),     intent(in)           :: directory
  type(String),     intent(in)           :: seedname
  real(dp),         intent(in), optional :: cutoff
  type(KpointGrid), intent(in), optional :: kpoint_grid
  real(dp),         intent(in), optional :: smearing
  
  type(IFile)       :: input_file
  type(QeInputFile) :: qe_input_file
  type(OFile)       :: output_file
  
  integer  :: ecutwfc_line
  integer  :: ecutrho_line
  real(dp) :: ecutwfc
  real(dp) :: ecutrho
  
  integer :: smearing_line
  
  type(String), allocatable :: line(:)
  
  integer :: i
  
  ! Read input file.
  input_file = IFile(seedname//'.in')
  qe_input_file = QeInputFile(input_file%lines())
  
  ! Update k-point grid.
  if (present(kpoint_grid)) then
    line = split_line(qe_input_file%k_points(2))
    if (size(line)==3) then
      qe_input_file%k_points(2) = str(kpoint_grid)
    elseif (size(line)==6) then
      qe_input_file%k_points(2) = join([str(kpoint_grid), line(4:6)])
    else
      call print_line(ERROR//': k_points card has an unexpected number of &
         &entries.')
      call quit()
    endif
  endif
  
  ! Update electronic cutoff.
  if (present(cutoff)) then
    ecutwfc_line = 0
    ecutrho_line = 0
    do i=1,size(qe_input_file%namelists)
      line = split_line(lower_case(qe_input_file%namelists(i)), delimiter='=')
      if (size(line)==2) then
        if (line(1)=='ecutwfc') then
          if (ecutwfc_line/=0) then
            call print_line(ERROR//': ecutwfc appears twice in input file.')
            call quit()
          endif
          ecutwfc_line = i
        elseif (line(1)=='ecutrho') then
          if (ecutrho_line/=0) then
            call print_line(ERROR//': ecutrho appears twice in input file.')
            call quit()
          endif
          ecutrho_line = i
        endif
      endif
    enddo
    
    if (ecutwfc_line==0) then
      call print_line(ERROR//': ecutwfc not present.')
      call err()
    endif
    
    line = split_line(qe_input_file%namelists(ecutwfc_line), delimiter='=')
    ecutwfc = dble(line(2))
    qe_input_file%namelists(ecutwfc_line) = &
       & 'ecutwfc='//cutoff*RYDBERG_PER_HARTREE//','
    
    if (ecutrho_line/=0) then
      line = split_line(qe_input_file%namelists(ecutrho_line), delimiter='=')
      ecutrho = dble(line(2))
      qe_input_file%namelists(ecutrho_line) = 'ecutrho='// &
         & ecutrho*cutoff/ecutwfc*RYDBERG_PER_HARTREE//','
    endif
  endif
  
  ! Update electronic smearing.
  if (present(smearing)) then
    smearing_line = 0
    do i=1,size(qe_input_file%namelists)
      line = split_line(qe_input_file%namelists(i), delimiter='=')
      if (line(1)=='degauss') then
        smearing_line = i
      endif
    enddo
    
    if (smearing_line==0) then
      call print_line(ERROR//': "degauss" not present in .in file.')
      call err()
    else
      qe_input_file%namelists(smearing_line) = 'degauss='// &
         & smearing*RYDBERG_PER_HARTREE//','
    endif
  endif
  
  ! Write output file.
  output_file = OFile(directory//'/'//seedname//'.in')
  call output_file%print_lines(qe_input_file)
end subroutine

! Find the maximum element-wise difference between two vectors.
impure elemental function maximum_difference(this,that) result(output)
  implicit none
  
  type(RealVector), intent(in) :: this
  type(RealVector), intent(in) :: that
  real(dp)                     :: output
  
  if (size(this)/=size(that)) then
    call print_line(CODE_ERROR//': Vectors of different lengths.')
    call err()
  endif
  
  output = maxval(abs(dble(this-that)))
end function
end module
