submodule (caesar_calculate_anharmonic_observables_module) caesar_calculate_anharmonic_observables_submodule
  use caesar_anharmonic_module
contains

module procedure calculate_anharmonic_observables_mode
  output%mode_name = 'calculate_anharmonic_observables'
  output%description = 'Calculates observables under the VSCF approximation. &
     &Should be run after calculate_potential.'
  output%keywords = [                                                         &
     & KeywordData( 'use_interpolated_potential', &
     &              'use_interpolated_potential specifies whether to use the &
     &interpolated potential from calculate_potential or to use the potential &
     &as calculated on the anharmonic q-point grid. If &
     &use_interpolated_potential is true then stress cannot be calculated.',  &
     &              default_value='false'),                                   &
     & KeywordData( 'split_q-points',                                         &
     &              'split_q-points specifies whether or not to split up VSCF &
     &subspaces by q-point. Doing so makes the calculation somewhat less &
     &accurate, but is required for reasonable calculation times if the &
     &system contains large subspaces spread over multiple q-points.',        &
     &              default_value='true'),                                    &
     & KeywordData( 'energy_convergence',                                     &
     &              'energy_convergence is the precision to which energies &
     &will be converged when constructing the VSCF ground state.' ),          &
     & KeywordData( 'no_converged_calculations_vscf',                         &
     &              'no_converged_calculations_vscf is the number of &
     &consecutive calculations which must be converged to within &
     &energy_convergence for the VSCF procedure to terminate.',               &
     &              default_value='5' ),                                      &
     & KeywordData( 'pre_pulay_iterations',                                   &
     &              'pre_pulay_iterations is the number of damped iterations &
     &which will be performed before the Pulay scheme is called. This must be &
     &at least 2.',                                                           &
     &              default_value='2' ),                                      &
     & KeywordData( 'pre_pulay_damping',                                      &
     &              "pre_pulay_damping is the damping factor of the pre-Pulay &
     &iterations. If potential at the start and end of iteration i are Vi and &
     &Vi' respectively, then the potential at the start of iteration j, Vj, &
     &is given by Vj = d*Vi+(1-d)*Vi'. Pre_pulay_damping must be between 0 &
     &and 1 inclusive.",                                                      &
     &              default_value='0.9' ),                                    &
     & KeywordData( 'max_pulay_iterations',                                   &
     &              'max_pulay_iterations is the maximum number of &
     &self-consistency iterations which will be passed into the Pulay &
     &scheme. This must be at least 2.',                                      &
     &              default_value='100' ),                                    &
     & KeywordData( 'iterative_mixing',                                       &
     &              'iterative_mixing controls how much the Pulay scheme is &
     &mixed with a damped iterative scheme. If iterative_mixing is 0, the &
     &iterations are entirely Pulay, and if it is 1 then the iterations are &
     &entirely damped iterative. iterative_mixing must be between 0 and 1 &
     &inclusive.',                                                            &
     &              default_value='0.1'),                                     &
     & KeywordData( 'iterative_damping',                                      &
     &              'iterative_damping is the damping factor of the iterative &
     &scheme which is mixed into the Pulay scheme. iterative_damping must be &
     &between 0 and 1 inclusive.',                                            &
     &              default_value='0.1'),                                     &
     & KeywordData( 'pulay_noise',                                            &
     &              "pulay_noise controls the random noise added to each &
     &Pulay iteration. If a Pulay iteration has coefficients xi and the &
     &noise-free next iteration would have coefficients xi', then the next &
     &iteration includes noise up to (xi'-xi)*pulay_noise. Noise is included &
     &to help the Pulay scheme break out of linearly-dependent subspaces.",   &
     &              default_value='1e-6'),                                    &
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin. min_temperature should be greater than zero, as the &
     &VSCF equations may not have a solution at zero temeperature.'),         &
     & KeywordData( 'max_temperature',                                        &
     &              'max_temperature is the maximum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.'),                                                     &
     & KeywordData( 'no_temperature_steps',                                   &
     &              'no_temperature_steps is the number of temperatures at &
     &which thermodynamic quantities are calculated.',                        &
     &              default_value='0'),                                       &
     & KeywordData( 'no_vscf_basis_states',                                   &
     &              'no_vscf_basis_states is the number of states along each &
     &mode in the basis used for the VSCF calculation.'),                     &
     & KeywordData( 'state_energy_cutoff',                                    &
     &              'state_energy_cutoff is the maximum value of <x|H|x> at &
     &which a basis state |x> is included in the calculation of the density &
     &matrix in the VSCF calculation. This energy is relative to the minimum &
     &value of <x|H|x> across all states |x>. state_energy_cutoff should be &
     &given in Hartree.'),                                                    &
     & KeywordData( 'min_frequency',                                          &
     &              'min_frequency is the frequency below which modes will be &
     &ignored when calculating thermodynamic quantities. min_frequency should &
     &be given in Hartree.',                                                  &
     &              default_value='1e-8'),                                    &
     & KeywordData( 'path',                                                   &
     &              'path is the path through fractional reciprocal space &
     &which will be mapped by the phonon dispersion curve. The path should be &
     &specified as labels and q-points, separated by commas. Breaks in the &
     &path should be marked with pipes, "|". The Gamma-point should be &
     &labelled G.',                                                           &
     &              default_value='G 0.0 0.0 0.0, R 0.5 0.5 0.5, &
     &M 0.0 0.5 0.5, G 0.0 0.0 0.0, X 0.0 0.0 0.5 | R 0.5 0.5 0.5, &
     &X 0.0 0.0 0.5, M 0.0 0.5 0.5'),                                         &
     & KeywordData( 'no_path_points',                                         &
     &              'no_path_points is the number of q-points sampled along &
     &the entire dispersion curve path.',                                     &
     &              default_value='1000'),                                    &
     & KeywordData( 'no_dos_samples',                                         &
     &              'no_dos_samples is the number of points in reciprocal &
     &space at which the normal modes are calculated when calculating the &
     &vibrational density of states.',                                        &
     &              default_value='10000')                                    ]
  output%main_subroutine => calculate_anharmonic_observables_subroutine
end procedure

module procedure calculate_anharmonic_observables_subroutine
  ! Inputs.
  type(RandomReal)      :: random_generator
  integer               :: random_generator_seed
  logical               :: use_interpolated_potential
  logical               :: split_qpoints
  real(dp)              :: energy_convergence
  integer               :: no_converged_calculations
  integer               :: pre_pulay_iterations
  real(dp)              :: pre_pulay_damping
  integer               :: max_pulay_iterations
  real(dp)              :: iterative_mixing
  real(dp)              :: iterative_damping
  real(dp)              :: pulay_noise
  real(dp)              :: min_temperature
  real(dp)              :: max_temperature
  integer               :: no_temperature_steps
  integer               :: no_vscf_basis_states
  real(dp)              :: state_energy_cutoff
  real(dp)              :: min_frequency
  real(dp), allocatable :: thermal_energies(:)
  type(String)          :: path
  integer               :: no_path_points
  integer               :: no_dos_samples
  
  ! Previous inputs.
  type(Dictionary)                   :: setup_harmonic_arguments
  type(String)                       :: seedname
  type(Dictionary)                   :: setup_anharmonic_arguments
  type(String)                       :: harmonic_path
  logical                            :: calculate_stress
  type(String)                       :: potential_representation
  type(StructureData)                :: harmonic_supercell
  type(QpointData),      allocatable :: harmonic_qpoints(:)
  type(DynamicalMatrix), allocatable :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     allocatable :: harmonic_complex_modes(:,:)
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(ComplexMode),        allocatable :: modes(:)
  type(QpointData),         allocatable :: qpoints(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  type(StructureData)                   :: anharmonic_supercell
  
  ! Convergence data.
  type(ConvergenceData) :: convergence_data
  
  ! Subspace variables.
  type(ComplexMode),      allocatable :: subspace_modes(:)
  type(QpointData),       allocatable :: subspace_qpoints(:)
  type(StressPrefactors), allocatable :: stress_prefactors(:)
  type(RealMatrix),       allocatable :: stress_prefactor
  class(StressData),      allocatable :: subspace_stress
  
  ! max(harmonic frequencies, frequency of max displacement).
  real(dp), allocatable :: starting_frequencies(:)
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Anharmonic stress.
  class(StressData), allocatable :: stress
  
  ! VSCHA basis, states, potential and stress.
  type(HarmonicBasis),    allocatable :: vscha_basis(:)
  type(VscfOutput),       allocatable :: vscha_output(:)
  type(HarmonicStates),   allocatable :: vscha_states(:)
  type(PotentialPointer), allocatable :: vscha_potentials(:)
  type(StressPointer),    allocatable :: vscha_stresses(:)
  real(dp),               allocatable :: vscha_frequencies(:,:)
  
  ! VSCF basis, states, potential and stress.
  type(SubspaceBasisPointer), allocatable :: vscf_basis(:)
  type(VscfOutput),           allocatable :: vscf_output(:)
  type(PotentialPointer),     allocatable :: vscf_potentials(:)
  type(BasisStatesPointer),   allocatable :: vscf_states(:)
  type(StressPointer),        allocatable :: vscf_stresses(:)
  
  ! Variables for interpolating VSCHA and VSCF.
  type(ComplexMode),           allocatable :: vscha_modes(:,:)
  type(ComplexMode),           allocatable :: vscf_modes(:,:)
  type(DynamicalMatrix),       allocatable :: dynamical_matrices(:)
  type(ComplexMode),           allocatable :: qpoint_modes(:)
  real(dp),                    allocatable :: qpoint_frequencies(:)
  type(CartesianHessian)                   :: hessian
  type(MinImages),             allocatable :: anharmonic_min_images(:,:)
  type(MinImages),             allocatable :: harmonic_min_images(:,:)
  type(StressDynamicalMatrix), allocatable :: stress_dynamical_matrices(:)
  type(StressHessian),         allocatable :: stress_hessian
  
  ! Dynamical matrix and hessian equivalents for the stress.
  
  ! Dispersion and density of states.
  type(PhononDispersion) :: phonon_dispersion
  type(PhononDos)        :: phonon_dos
  
  ! Corrections for double-counting of energies and stresses.
  real(dp)         :: energy_correction
  type(RealMatrix) :: stress_correction
  
  ! Thermodynamic data.
  type(ThermodynamicData), allocatable :: interpolated_vscha_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscha1_thermodynamics(:,:)
  type(ThermodynamicData), allocatable :: vscha2_thermodynamics(:,:)
  type(ThermodynamicData), allocatable :: vscf_thermodynamics(:,:)
  type(ThermodynamicData), allocatable :: interpolated_vscf_thermodynamics(:)
  
  ! Input files.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: potential_file
  type(IFile)  :: stress_file
  type(IFile)  :: harmonic_supercell_file
  type(IFile)  :: harmonic_qpoints_file
  type(String) :: qpoint_dir
  type(IFile)  :: harmonic_dynamical_matrices_file
  type(IFile)  :: harmonic_complex_modes_file
  
  ! Output directories.
  type(String) :: observables_dir
  type(String) :: vscha_dir
  type(String) :: vscha_temperature_dir
  type(String) :: vscf_dir
  type(String) :: vscf_temperature_dir
  
  ! Output files.
  type(OFile) :: vscha_logfile
  type(OFile) :: vscf_logfile
  type(OFile) :: temperature_file
  type(OFile) :: convergence_file
  
  ! Effective harmonic output files.
  type(OFile) :: castep_phonon_file
  type(OFile) :: qe_force_constants_file
  
  ! Thermodynamic data output files.
  type(OFile) :: vscha1_thermodynamic_file
  type(OFile) :: vscha2_thermodynamic_file
  type(OFile) :: vscf_thermodynamic_file
  type(OFile) :: interpolated_vscha_thermodynamic_file
  type(OFile) :: interpolated_vscf_thermodynamic_file
  
  ! Working variables.
  type(String) :: thermodynamic_header
  logical      :: can_be_interpolated
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in arguments from user.
  ! --------------------------------------------------
  if (arguments%is_set('random_seed')) then
    random_generator = RandomReal(int(arguments%value('random_seed')))
  else
    random_generator = RandomReal()
  endif
  random_generator_seed = random_generator%get_seed()
  use_interpolated_potential = lgcl(                 &
     & arguments%value('use_interpolated_potential') )
  split_qpoints = lgcl(arguments%value('split_q-points'))
  energy_convergence = dble(arguments%value('energy_convergence'))
  no_converged_calculations = int(                       &
     & arguments%value('no_converged_calculations_vscf') )
  pre_pulay_iterations = int(arguments%value('pre_pulay_iterations'))
  pre_pulay_damping = dble(arguments%value('pre_pulay_damping'))
  max_pulay_iterations = int(arguments%value('max_pulay_iterations'))
  iterative_mixing = dble(arguments%value('iterative_mixing'))
  iterative_damping = dble(arguments%value('iterative_damping'))
  pulay_noise = dble(arguments%value('pulay_noise'))
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path = arguments%value('path')
  no_path_points = int(arguments%value('no_path_points'))
  no_dos_samples = int(arguments%value('no_dos_samples'))
  no_vscf_basis_states = int(arguments%value('no_vscf_basis_states'))
  state_energy_cutoff = dble(arguments%value('state_energy_cutoff'))
  
  ! Check inputs.
  if (min_temperature<0) then
    call print_line(ERROR//': min_temperature must not be less than 0K.')
    call quit()
  elseif (max_temperature<min_temperature) then
    call print_line(ERROR//': max_temperature must not be less than &
       &min_temperature.')
    call quit()
  elseif (min_frequency<0) then
    call print_line(ERROR//': min_frequency must not be less than 0 Hartree.')
    call quit()
  elseif (no_vscf_basis_states<1) then
    call print_line(ERROR//': no_vscf_basis_states must be at least 1.')
    call quit()
  elseif (pre_pulay_damping<0 .or. pre_pulay_damping>1) then
    call print_line(ERROR//': pre_pulay_damping must be between 0 and 1.')
    call quit()
  elseif (pre_pulay_iterations<2) then
    call print_line(ERROR//': pre_pulay_iterations must be at least 2.')
    call quit()
  elseif (max_pulay_iterations<2) then
    call print_line(ERROR//': max_pulay_iterations must be at least 2.')
    call quit()
  elseif (iterative_mixing<0 .or. iterative_mixing>1) then
    call print_line(ERROR//': iterative_mixing must be between 0 and 1.')
    call quit()
  elseif (iterative_damping<0 .or. iterative_damping>1) then
    call print_line(ERROR//': iterative_damping must be between 0 and 1.')
    call quit()
  endif
  
  if (min_temperature<=0) then
    call print_line(WARNING//': Calculation may not converge if &
       &min_temperature is exactly 0K.')
  endif
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = Dictionary(setup_harmonic_mode())
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = Dictionary(setup_anharmonic_mode())
  call setup_anharmonic_arguments%read_file('setup_anharmonic.used_settings')
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  calculate_stress = lgcl(setup_anharmonic_arguments%value('calculate_stress'))
  potential_representation = setup_anharmonic_arguments%value( &
                                  & 'potential_representation' )
  
  ! TODO: remove this constraint by interpolating stress correctly.
  if (use_interpolated_potential) then
    calculate_stress = .false.
  endif
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  
  ! Read in anharmonic data.
  if (use_interpolated_potential) then
    anharmonic_data_file = IFile('interpolated_anharmonic_data.dat')
    potential_file = IFile('interpolated_potential.dat')
  else
    anharmonic_data_file = IFile('anharmonic_data.dat')
    potential_file = IFile('potential.dat')
  endif
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  potential = PotentialPointer(potential_file%lines())
  
  modes = anharmonic_data%complex_modes
  qpoints = anharmonic_data%qpoints
  subspaces = anharmonic_data%degenerate_subspaces
  anharmonic_supercell = anharmonic_data%anharmonic_supercell
  
  ! Read in anharmonic stress.
  if (calculate_stress) then
    stress_file = IFile('stress.dat')
    stress = StressPointer(stress_file%lines())
  endif
  
  ! Read in harmonic q-points and dynamical matrices.
  harmonic_supercell_file = IFile(harmonic_path//'/large_supercell.dat')
  harmonic_supercell = StructureData(harmonic_supercell_file%lines())
  harmonic_qpoints_file = IFile(harmonic_path//'/qpoints.dat')
  harmonic_qpoints = QpointData(harmonic_qpoints_file%sections())
  allocate( harmonic_dynamical_matrices(size(harmonic_qpoints)),           &
          & harmonic_complex_modes( anharmonic_data%structure%no_modes,    &
          &                         size(harmonic_qpoints)              ), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_qpoints)
    qpoint_dir = harmonic_path// &
               & '/qpoint_'//left_pad(i,str(size(harmonic_qpoints)))
    harmonic_dynamical_matrices_file = IFile( &
        & qpoint_dir//'/dynamical_matrix.dat' )
    harmonic_dynamical_matrices(i) = DynamicalMatrix( &
           & harmonic_dynamical_matrices_file%lines() )
    
    harmonic_complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    harmonic_complex_modes(:,i) = ComplexMode(  &
       & harmonic_complex_modes_file%sections() )
  enddo
  
  ! --------------------------------------------------
  ! Generate objects which don't depend on the potential.
  ! --------------------------------------------------
  ! Generate thermal energies.
  ! thermal_energies(1)                    = min_temperature * kB.
  ! thermal_energies(no_temperature_steps) = max_temperature * kB.
  if (no_temperature_steps==1) then
    thermal_energies = [KB_IN_AU*min_temperature]
  else
    allocate( thermal_energies(no_temperature_steps), &
            & stat=ialloc); call err(ialloc)
    do i=1,no_temperature_steps
      thermal_energies(i) = KB_IN_AU                                   &
                        & * ( min_temperature*(no_temperature_steps-i) &
                        &   + max_temperature*(i-1) )                  &
                        & / (no_temperature_steps-1.0_dp)
    enddo
  endif
  
  ! Construct minimum image data.
  harmonic_min_images = calculate_min_images(harmonic_supercell)
  anharmonic_min_images = calculate_min_images(anharmonic_supercell)
  
  ! Generate stress prefactors, which are geometric factors quantifying
  !    the direction of each pair of modes w/r/t the lattice vectors.
  if (calculate_stress) then
    stress_prefactors = [( StressPrefactors(subspaces(i),modes), &
                         & i=1,                                  &
                         & size(subspaces)                       )]
  endif
  
  ! Calculate starting frequencies.
  ! If the subspace frequencies are large enough, then the starting frequency
  !    is simply the harmonic frequency.
  ! The starting frequencies must be at least equal to a number of values
  !    to ensure that convergence isn't falsely reached early
  !    because the frequency changes with each step are too small.
  starting_frequencies = [max(                                    &
     & subspaces%frequency,                                       &
     & maxval([ anharmonic_data%frequency_of_max_displacement,    &
     &          min_frequency,                                    &
     &          energy_convergence                             ]) )]
  
  ! Construct convergence data.
  convergence_data = ConvergenceData( pre_pulay_iterations,     &
                                    & pre_pulay_damping,        &
                                    & max_pulay_iterations,     &
                                    & iterative_mixing,         &
                                    & iterative_damping,        &
                                    & pulay_noise,              &
                                    & energy_convergence,       &
                                    & no_converged_calculations )
  
  ! --------------------------------------------------
  ! Generate directories.
  ! --------------------------------------------------
  ! Make anharmonic_observables directory, and VSCHA and VSCF subdirectories.
  observables_dir = 'anharmonic_observables'
  call mkdir(observables_dir)
  
  vscha_dir = observables_dir//'/vscha'
  call mkdir(vscha_dir)
  vscha_logfile = OFile(vscha_dir//'/logfile.dat')
  
  vscf_dir = observables_dir//'/vscf'
  call mkdir(vscf_dir)
  vscf_logfile = OFile(vscf_dir//'/logfile.dat')
      
  
  ! Write temperature file in the VSCHA and VSCF directories.
  temperature_file = OFile(vscha_dir//'/temperatures.dat')
  call temperature_file%print_line('Thermal energies, KbT, (Ha):')
  call temperature_file%print_lines(thermal_energies)
  
  temperature_file = OFile(vscf_dir//'/temperatures.dat')
  call temperature_file%print_line('Thermal energies, KbT, (Ha):')
  call temperature_file%print_lines(thermal_energies)
  
  ! --------------------------------------------------
  ! Run calculations at each temperature.
  ! --------------------------------------------------
  allocate( vscha_frequencies( size(subspaces),                        &
          &                    size(thermal_energies) ),               &
          & vscf_basis(size(subspaces)),                               &
          & vscha_modes( anharmonic_data%structure%no_modes,           &
          &              size(qpoints)                       ),        &
          & dynamical_matrices(size(qpoints)),                         &
          & vscf_modes( anharmonic_data%structure%no_modes,            &
          &             size(qpoints)                       ),         &
          & vscf_thermodynamics( size(subspaces),                      &
          &                      size(thermal_energies) ),             &
          & vscha1_thermodynamics( size(subspaces),                    &
          &                        size(thermal_energies) ),           &
          & vscha2_thermodynamics( size(subspaces),                    &
          &                        size(thermal_energies) ),           &
          & interpolated_vscha_thermodynamics(size(thermal_energies)), &
          & interpolated_vscf_thermodynamics(size(thermal_energies)),  &
          & stat=ialloc); call err(ialloc)
  do i=1,size(thermal_energies)
    call print_line('')
    call print_line(repeat('=',50))
    call print_line('Thermal energy '//i//' of '//size(thermal_energies)// &
       & ': KT = '//thermal_energies(i)//' Ha')
    call print_line(repeat('=',50))
    
    ! --------------------------------------------------
    ! Run VSCF using VSCHA basis.
    ! --------------------------------------------------
    
    ! Make VSCHA temperature directory.
    vscha_temperature_dir = vscha_dir//'/temperature_'// &
       & left_pad(i,str(size(thermal_energies)))
    call mkdir(vscha_temperature_dir)
    
    if (i==1) then
      vscha_basis = HarmonicBasis( subspaces%id,                &
                                 & starting_frequencies,        &
                                 & anharmonic_supercell%sc_size )
    else
      vscha_basis = HarmonicBasis( subspaces%id,                &
                                 & vscha_frequencies(:,i-1),    &
                                 & anharmonic_supercell%sc_size )
    endif
    
    call print_line('Running VSCHA.')
    ! N.B. for the first run, vscha_output is not allocated,
    !    so starting_configuration will not be present() inside run_vscf.
    ! For later runs, starting_configuration will be the output of the
    !    previous run.
    convergence_file = OFile( vscha_temperature_dir// &
                            & '/convergence.dat'      )
    vscha_output = run_vscf( potential,                                &
                           & stress,                                   &
                           & subspaces,                                &
                           & vscha_basis,                              &
                           & thermal_energies(i),                      &
                           & state_energy_cutoff,                      &
                           & starting_frequencies,                     &
                           & convergence_data,                         &
                           & anharmonic_data,                          &
                           & random_generator,                         &
                           & starting_configuration = vscha_output,    &
                           & convergence_file       = convergence_file )
    vscha_potentials = [(vscha_output(j)%potential, j=1, size(vscha_output))]
    vscha_states = [( HarmonicStates(vscha_output(i)%states), &
                    & i=1,                                    &
                    & size(vscha_output)                      )]
    vscha_frequencies(:,i) = vscha_states%frequency
    if (any(vscha_states%frequency<=energy_convergence)) then
      call print_line(WARNING//': At least one VSCHA frequency is below &
         &energy_convergence. Consider lowering energy_convergence.')
      call print_line('Minimum frequency: '//minval(vscha_states%frequency))
    endif
    
    if (calculate_stress) then
      vscha_stresses = [(vscha_output(j)%stress, j=1, size(vscha_output))]
    endif
    call print_line('')
    
    !! Calculate thermodynamic quantities with VSCHA basis.
    call print_line('Calculating thermodynamic observables under VSCHA.')
    do j=1,size(subspaces)
      if (calculate_stress) then
        subspace_stress = vscha_stresses(j)
        stress_prefactor = stress_prefactors(j)%average_prefactor()
      endif
      
      vscha1_thermodynamics(j,i) = harmonic_observables(       &
         &    thermal_energy   = thermal_energies(i),          &
         &    stress           = subspace_stress,              &
         &    stress_prefactor = stress_prefactor,             &
         &    frequency        = vscha_frequencies(j,i),       &
         &    num_dimensions   = size(subspaces(j)),           &
         &    supercell_size   = anharmonic_supercell%sc_size, &
         &    anharmonic_data  = anharmonic_data               ) &
         & / anharmonic_supercell%sc_size
      
      vscha2_thermodynamics(j,i) = effective_harmonic_observables( &
           &    thermal_energy   = thermal_energies(i),            &
           &    potential        = vscha_potentials(j),            &
           &    stress           = subspace_stress,                &
           &    stress_prefactor = stress_prefactor,               &
           &    frequency        = vscha_frequencies(j,i),         &
           &    num_dimensions   = size(subspaces(j)),             &
           &    supercell_size   = anharmonic_supercell%sc_size,   &
           &    anharmonic_data  = anharmonic_data               ) &
           & / anharmonic_supercell%sc_size
    enddo
    
    ! --------------------------------------------------
    ! Calculated interpolated VSCHA results.
    ! --------------------------------------------------
    
    call print_line('Interpolating VSCHA results.')
    
    ! Assemble the dynamical matrices at each q-point from the complex modes
    !    and the effective frequencies.
    do j=1,size(qpoints)
      qpoint_modes = modes(filter(modes%qpoint_id==qpoints(j)%id))
      qpoint_frequencies = [(                                      &
         & vscha_frequencies(                                      &
         &    first(subspaces%id==qpoint_modes(k)%subspace_id),    &
         &    i                                                 ), &
         & k=1,                                                    &
         & size(qpoint_modes)                                      )]
      qpoint_frequencies = max(qpoint_frequencies, min_frequency)
      if (size(qpoint_modes)==anharmonic_data%structure%no_modes) then
        vscha_modes(:,j) = qpoint_modes
        vscha_modes(:,j)%frequency = qpoint_frequencies
      else
        vscha_modes(:3,j) = generate_translational_modes( &
                             & anharmonic_data%structure, &
                             & qpoints                    )
        vscha_modes(4:,j) = qpoint_modes
        vscha_modes(4:,j)%frequency = qpoint_frequencies
      endif
      dynamical_matrices(j) = DynamicalMatrix(qpoint_modes,qpoint_frequencies)
    enddo
    
    ! Construct Hessian from dynamical matrices.
    hessian = reconstruct_hessian( anharmonic_supercell, &
                                 & qpoints,              &
                                 & dynamical_matrices,   &
                                 & vscha_logfile         )
    
    ! Write out Castep .phonon file and QE .fc file for the VSCHA modes.
    castep_phonon_file = OFile( vscha_temperature_dir                 &
                      & //'/'// make_castep_phonon_filename(seedname) )
    call write_castep_phonon_file( castep_phonon_file,       &
                                 & vscha_modes,              &
                                 & qpoints,                  &
                                 & anharmonic_data%structure )
    
    qe_force_constants_file = OFile(                        &
       &         vscha_temperature_dir                      &
       & //'/'// make_qe_force_constants_filename(seedname) )
    call write_qe_force_constants_file( qe_force_constants_file,   &
                                      & hessian,                   &
                                      & anharmonic_data%structure, &
                                      & anharmonic_supercell       )
    
    ! Interpolate hessian.
    hessian = calculate_interpolated_hessian( anharmonic_supercell,        &
                                            & qpoints,                     &
                                            & dynamical_matrices,          &
                                            & anharmonic_min_images,       &
                                            & harmonic_supercell,          &
                                            & harmonic_qpoints,            &
                                            & harmonic_dynamical_matrices, &
                                            & harmonic_min_images          )
      
    ! Calculate stress.
    if (calculate_stress) then
      stress_dynamical_matrices = stress%calculate_dynamical_matrices( &
                                                & qpoints,             &
                                                & thermal_energies(i), &
                                                & subspaces,           &
                                                & vscha_basis,         &
                                                & vscha_states,        &
                                                & anharmonic_data      )
      
      stress_hessian = reconstruct_stress_hessian( &
                      & anharmonic_supercell,      &
                      & qpoints,                   &
                      & stress_dynamical_matrices, &
                      & vscha_logfile              )
    endif
    
    ! Generate self-consistent phonon dispersion curve by interpolating between
    !    calculated q-points using Fourier interpolation.
    phonon_dispersion = PhononDispersion(           &
       & supercell         = harmonic_supercell,    &
       & min_images        = harmonic_min_images,   &
       & hessian           = hessian,               &
       & stress_hessian    = stress_hessian,        &
       & stress_supercell  = anharmonic_supercell,  &
       & stress_min_images = anharmonic_min_images, &
       & path_string       = path,                  &
       & no_path_points    = no_path_points,        &
       & logfile           = vscha_logfile          )
    call phonon_dispersion%write_files( vscha_temperature_dir,    &
                                      & seedname,                 &
                                      & anharmonic_data%structure )
    
    ! Generate self-consistent phonon density of states,
    !    interpolating as above.
    ! Re-seed the random generator each time, so that the results at different
    !    temperatures can be compared free from random noise.
    random_generator = RandomReal(random_generator_seed)
    phonon_dos = PhononDos( supercell         = harmonic_supercell,    &
                          & min_images        = harmonic_min_images,   &
                          & hessian           = hessian,               &
                          & stress_hessian    = stress_hessian,        &
                          & stress_supercell  = anharmonic_supercell,  &
                          & stress_min_images = anharmonic_min_images, &
                          & thermal_energies  = [thermal_energies(i)], &
                          & min_frequency     = min_frequency,         &
                          & no_dos_samples    = no_dos_samples,        &
                          & logfile           = vscha_logfile,         &
                          & random_generator  = random_generator       )
    call phonon_dos%write_files(vscha_temperature_dir)
    
    interpolated_vscha_thermodynamics(i) = phonon_dos%thermodynamic_data(1)
    
    ! --------------------------------------------------
    ! Run VSCF using VCI basis.
    ! --------------------------------------------------
    
    ! Make VSCF temperature directory.
    vscf_temperature_dir = vscf_dir//'/temperature_'// &
       & left_pad(i,str(size(thermal_energies)))
    call mkdir(vscf_temperature_dir)
    
    ! Generate basis of states.
    call print_line('')
    call print_line('Generating basis for VSCF.')
    do j=1,size(subspaces)
      subspace_modes = subspaces(j)%modes(modes)
      subspace_qpoints = subspaces(j)%qpoints(modes,qpoints)
      subspace_qpoints = subspace_qpoints(set(subspace_qpoints%id))
      call print_line( 'Generating basis in subspace ' // &
                     & subspaces(j)%id                 // &
                     & ', containing '                 // &
                     & size(subspace_modes)            // &
                     & ' modes across '                // &
                     & size(subspace_qpoints)          // &
                     & ' q-points.'                       )
      if (split_qpoints) then
        vscf_basis(j) = SubspaceBasisPointer(SplitQpointsBasis( &
              & subspaces(j),                                   &
              & vscha_frequencies(j,i),                         &
              & modes,                                          &
              & qpoints,                                        &
              & anharmonic_supercell,                           &
              & no_vscf_basis_states-1,                         &
              & anharmonic_data%potential_expansion_order,      &
              & anharmonic_data%structure%symmetries            ))
      else
        vscf_basis(j) = SubspaceBasisPointer(FullSubspaceBasis( &
              & subspaces(j),                                   &
              & vscha_frequencies(j,i),                         &
              & modes,                                          &
              & qpoints,                                        &
              & anharmonic_supercell,                           &
              & no_vscf_basis_states-1,                         &
              & anharmonic_data%potential_expansion_order,      &
              & anharmonic_data%structure%symmetries            ))
      endif
    enddo
    call print_line('')
    
    ! Run VSCF to generate single-subspace potentials and states.
    ! N.B. for the first run, vscf_output is not allocated,
    !    so starting_configuration will not be present() inside run_vscf.
    ! For later runs, starting_configuration will be the output of the
    !    previous run.
    call print_line('Running VSCF.')
    convergence_file = OFile(vscf_temperature_dir//'/convergence.dat')
    vscf_output = run_vscf( potential,                                &
                          & stress,                                   &
                          & subspaces,                                &
                          & vscf_basis,                               &
                          & thermal_energies(i),                      &
                          & state_energy_cutoff,                      &
                          & vscha_frequencies(:,i),                   &
                          & convergence_data,                         &
                          & anharmonic_data,                          &
                          & random_generator,                         &
                          & starting_configuration = vscf_output,     &
                          & convergence_file       = convergence_file )
    
    vscf_potentials = [(vscf_output(j)%potential, j=1, size(vscf_output))]
    vscf_states = [(vscf_output(j)%states, j=1, size(vscf_output))]
    if (calculate_stress) then
      vscf_stresses = [(vscf_output(j)%stress, j=1, size(vscf_output))]
    endif
    call print_line('')
    
    ! Calculate thermodynamic quantities with VCI basis.
    call print_line('Calculating thermodynamic observables under VSCF.')
    if (calculate_stress) then
      vscf_thermodynamics(:,i) = vscf_basis%thermodynamic_data( &
                             &           thermal_energies(i),   &
                             &           vscf_states,           &
                             &           subspaces,             &
                             &           vscf_potentials,       &
                             &           vscf_stresses,         &
                             &           stress_prefactors,     &
                             &           anharmonic_data      ) &
                             & / anharmonic_supercell%sc_size
    else
      vscf_thermodynamics(:,i) = vscf_basis%thermodynamic_data(              &
                             &           thermal_energies(i),                &
                             &           vscf_states,                        &
                             &           subspaces,                          &
                             &           vscf_potentials,                    &
                             &           anharmonic_data = anharmonic_data ) &
                             & / anharmonic_supercell%sc_size
    endif
    
    ! --------------------------------------------------
    ! Calculated interpolated VSCF results.
    ! --------------------------------------------------
    if (.not. calculate_stress) then
      can_be_interpolated = potential%can_be_interpolated()
    else
      can_be_interpolated = potential%can_be_interpolated() &
                    & .and. stress%can_be_interpolated()
    endif
    
    if (can_be_interpolated) then
      call print_line('Interpolating VSCF results.')
      dynamical_matrices = potential%calculate_dynamical_matrices( &
                                            & qpoints,             &
                                            & thermal_energies(i), &
                                            & subspaces,           &
                                            & vscf_basis,          &
                                            & vscf_states,         &
                                            & anharmonic_data      )
      do j=1,size(dynamical_matrices)
        call dynamical_matrices(j)%check( anharmonic_data%structure, &
                                        & vscf_logfile               )
      enddo
      
      ! Construct effective modes from dynamical matrices.
      do j=1,size(qpoints)
        vscf_modes(:,j) = ComplexMode(                &
           & dynamical_matrices(j),                   &
           & anharmonic_data%structure,               &
           & modes_real=qpoints(j)%is_paired_qpoint() )
      enddo
      
      ! Construct Hessian from dynamical matrices.
      hessian = reconstruct_hessian( anharmonic_supercell, &
                                   & qpoints,              &
                                   & dynamical_matrices,   &
                                   & vscf_logfile          )
      
      ! Write out Castep .phonon file and QE .fc file for the VSCHA modes.
      castep_phonon_file = OFile( vscf_temperature_dir                 &
                        & //'/'// make_castep_phonon_filename(seedname) )
      call write_castep_phonon_file( castep_phonon_file,       &
                                   & vscf_modes,               &
                                   & qpoints,                  &
                                   & anharmonic_data%structure )
      
      qe_force_constants_file = OFile(                        &
         &         vscf_temperature_dir                       &
         & //'/'// make_qe_force_constants_filename(seedname) )
      call write_qe_force_constants_file( qe_force_constants_file,   &
                                        & hessian,                   &
                                        & anharmonic_data%structure, &
                                        & anharmonic_supercell       )
      
      ! Interpolate hessian.
      hessian = calculate_interpolated_hessian( anharmonic_supercell,        &
                                              & qpoints,                     &
                                              & dynamical_matrices,          &
                                              & anharmonic_min_images,       &
                                              & harmonic_supercell,          &
                                              & harmonic_qpoints,            &
                                              & harmonic_dynamical_matrices, &
                                              & harmonic_min_images          )
      
      ! Calculate stress.
      if (calculate_stress) then
        stress_dynamical_matrices = stress%calculate_dynamical_matrices( &
                                                  & qpoints,             &
                                                  & thermal_energies(i), &
                                                  & subspaces,           &
                                                  & vscf_basis,          &
                                                  & vscf_states,         &
                                                  & anharmonic_data      )
        stress_hessian = reconstruct_stress_hessian( &
                        & anharmonic_supercell,      &
                        & qpoints,                   &
                        & stress_dynamical_matrices, &
                        & vscf_logfile               )
      endif
      
      ! Generate self-consistent phonon dispersion curve by interpolating
      !    between calculated q-points using Fourier interpolation.
      phonon_dispersion = PhononDispersion(           &
         & supercell         = harmonic_supercell,    &
         & min_images        = harmonic_min_images,   &
         & hessian           = hessian,               &
         & stress_hessian    = stress_hessian,        &
         & stress_supercell  = anharmonic_supercell,  &
         & stress_min_images = anharmonic_min_images, &
         & path_string       = path,                  &
         & no_path_points    = no_path_points,        &
         & logfile           = vscf_logfile           )
      call phonon_dispersion%write_files( vscf_temperature_dir,     &
                                        & seedname,                 &
                                        & anharmonic_data%structure )
      
      ! Generate self-consistent phonon density of states,
      !    interpolating as above.
      ! Re-seed the random generator each time, so that the results at
      !    different temperatures can be compared free from random noise.
      random_generator = RandomReal(random_generator_seed)
      phonon_dos = PhononDos( supercell         = harmonic_supercell,    &
                            & min_images        = harmonic_min_images,   &
                            & hessian           = hessian,               &
                            & stress_hessian    = stress_hessian,        &
                            & stress_supercell  = anharmonic_supercell,  &
                            & stress_min_images = anharmonic_min_images, &
                            & thermal_energies  = [thermal_energies(i)], &
                            & min_frequency     = min_frequency,         &
                            & no_dos_samples    = no_dos_samples,        &
                            & logfile           = vscf_logfile,          &
                            & random_generator  = random_generator       )
      call phonon_dos%write_files(vscf_temperature_dir)
      
      ! Record thermodynamic data.
      interpolated_vscf_thermodynamics(i) = phonon_dos%thermodynamic_data(1)
      
      energy_correction = potential%energy_correction( subspaces,           &
                                                     & vscf_basis,          &
                                                     & vscf_states,         &
                                                     & anharmonic_data      )
      call interpolated_vscf_thermodynamics(i)%add_energy( &
          & energy_correction/anharmonic_supercell%sc_size )
      
      if (calculate_stress) then
        stress_correction = stress%stress_correction( subspaces,           &
                                                    & vscf_basis,          &
                                                    & vscf_states,         &
                                                    & anharmonic_data      )
        call interpolated_vscf_thermodynamics(i)%add_stress( &
            & stress_correction/anharmonic_supercell%sc_size )
      endif
    else
      if (.not. potential%can_be_interpolated()) then
        call print_line('Skipping VSCF interpolation as the chosen potential &
           &representation cannot be interpolated.')
      else
        call print_line('Skipping VSCF interpolation as the chosen stress &
           &representation cannot be interpolated.')
      endif
    endif
  enddo
  
  ! --------------------------------------------------
  ! Write out thermodynamic results.
  ! --------------------------------------------------
  
  thermodynamic_header =  &
     &'kB * temperature, T (Hartree) | &
     &Vibrational Energy per cell, U=<E>, (Hartree) | &
     &Vibrational Free Energy per cell, F=U-TS, (Hartree) | &
     &Vibrational Shannon Entropy per cell, S/k_B, (arb. units)'
  if (calculate_stress) then
    thermodynamic_header = thermodynamic_header // ' | &
       &Vibrational stress tensor (Hartree per bohr^3) (nine values) | &
       &Volume per cell, V (bohr^3 per cell) | &
       &Vibrational Enthalpy per cell, H=U+pV, (Hartree) | &
       &Vibrational Gibbs Free Energy per cell, G=U-TS+pV, (Hartree)'
  endif
  
  ! VSCHA without Fourier interpolation.
  vscha1_thermodynamic_file = OFile(                         &
     & observables_dir//'/vscha_thermodynamic_variables.dat' )
  call vscha1_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscha1_thermodynamic_file%print_line(sum(vscha1_thermodynamics(:,i)))
  enddo
  
  ! VSCHA without Fourier interpolation, but with <V> taken using the VSCF
  !    potential rather than the effective harmonic potential.
  vscha2_thermodynamic_file = OFile(                              &
     & observables_dir//'/vscha_vscf_thermodynamic_variables.dat' )
  call vscha2_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscha2_thermodynamic_file%print_line(sum(vscha2_thermodynamics(:,i)))
  enddo
  
  ! VSCF without Fourier interpolation.
  vscf_thermodynamic_file = OFile(                          &
     & observables_dir//'/vscf_thermodynamic_variables.dat' )
  call vscf_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscf_thermodynamic_file%print_line(sum(vscf_thermodynamics(:,i)))
  enddo
  
  ! VSCHA with Fourier interpolation.
  interpolated_vscha_thermodynamic_file = OFile(                          &
     & observables_dir//'/interpolated_vscha_thermodynamic_variables.dat' )
  call interpolated_vscha_thermodynamic_file%print_line(thermodynamic_header)
  call interpolated_vscha_thermodynamic_file%print_lines( &
                      & interpolated_vscha_thermodynamics )
  
  ! VSCF with Fourier interpolation.
  interpolated_vscf_thermodynamic_file = OFile(                          &
     & observables_dir//'/interpolated_vscf_thermodynamic_variables.dat' )
  call interpolated_vscf_thermodynamic_file%print_line(thermodynamic_header)
  call interpolated_vscf_thermodynamic_file%print_lines( &
                      & interpolated_vscf_thermodynamics )
end procedure
end submodule
