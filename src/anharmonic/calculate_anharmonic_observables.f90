! ======================================================================
! Calculates observables under the VSCF approximation.
! ======================================================================
! Should be run after calculate_potential.
module calculate_anharmonic_observables_module
  use common_module 
  use states_module
  use anharmonic_common_module
  use potentials_module
  
  use generate_subspace_potentials_module
  use vscf_module
  use stress_prefactors_module
  
  ! TODO: remove this.
  use interpolation_module
  implicit none
  
  private
  
  public :: startup_calculate_anharmonic_observables
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_calculate_anharmonic_observables()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'calculate_anharmonic_observables'
  mode%description = 'Calculates observables under the VSCF approximation. &
     &Should be run after calculate_potential.'
  mode%keywords = [                                                           &
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
     & KeywordData( 'max_pulay_iterations',                                   &
     &              'max_pulay_iterations is the maximum number of &
     &self-consistency iterations which will be passed into the Pulay &
     &scheme. This must be at least 2.',                                      &
     &              default_value='20' ),                                     &
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
     & KeywordData( 'min_temperature',                                        &
     &              'min_temperature is the minimum temperature at which &
     &thermodynamic quantities are calculated. min_temperature should be &
     &given in Kelvin.'),                                                     &
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
  mode%main_subroutine => calculate_anharmonic_observables_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! The main subroutine.
! ----------------------------------------------------------------------
subroutine calculate_anharmonic_observables_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Inputs.
  type(RandomReal)      :: random_generator
  integer               :: random_generator_seed
  logical               :: split_qpoints
  real(dp)              :: energy_convergence
  integer               :: no_converged_calculations_vscf
  integer               :: max_pulay_iterations
  integer               :: pre_pulay_iterations
  real(dp)              :: pre_pulay_damping
  real(dp)              :: min_temperature
  real(dp)              :: max_temperature
  integer               :: no_temperature_steps
  integer               :: no_vscf_basis_states
  real(dp)              :: min_frequency
  real(dp), allocatable :: thermal_energies(:)
  type(String)          :: path
  integer               :: no_path_points
  integer               :: no_dos_samples
  
  ! Previous inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(Dictionary) :: setup_anharmonic_arguments
  logical          :: calculate_stress
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(ComplexMode),        allocatable :: modes(:)
  type(QpointData),         allocatable :: qpoints(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  type(StructureData)                   :: supercell
  
  ! Subspace variables.
  type(ComplexMode),      allocatable :: subspace_modes(:)
  type(QpointData),       allocatable :: subspace_qpoints(:)
  type(StressPrefactors), allocatable :: stress_prefactors(:)
  type(RealMatrix),       allocatable :: stress_prefactor
  type(StressPointer),    allocatable :: subspace_stress
  
  ! max(harmonic frequencies, frequency of max displacement).
  real(dp), allocatable :: starting_frequencies(:)
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Anharmonic stress.
  type(StressPointer), allocatable :: stress
  
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
  type(PotentialPointer),     allocatable :: subspace_potentials(:)
  type(BasisStatesPointer),   allocatable :: subspace_states(:)
  type(StressPointer),        allocatable :: subspace_stresses(:)
  
  ! Interpolated VSCHA variables; modes and dynamical matrices at each q-point.
  type(ComplexMode),     allocatable :: vscha_modes(:,:)
  type(DynamicalMatrix), allocatable :: dynamical_matrices(:)
  type(ComplexMode),     allocatable :: qpoint_modes(:)
  real(dp),              allocatable :: qpoint_frequencies(:)
  
  ! Hessian matrix and minimum image data.
  type(CartesianHessian)       :: hessian
  type(MinImages), allocatable :: min_images(:,:)
  
  ! Dispersion and density of states.
  type(PhononDispersion) :: phonon_dispersion
  type(PhononDos)        :: phonon_dos
  
  ! Thermodynamic data.
  type(ThermodynamicData), allocatable :: vscha_thermodynamics(:)
  type(ThermodynamicData), allocatable :: vscha1_thermodynamics(:,:)
  type(ThermodynamicData), allocatable :: vscha2_thermodynamics(:,:)
  type(ThermodynamicData), allocatable :: vscf_thermodynamics(:,:)
  
  type(SubspaceWavefunctionsPointer), allocatable :: subspace_wavefunctions(:)
  
  type(String) :: thermodynamic_header
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: potential_file
  type(IFile)  :: stress_file
  type(String) :: output_dir
  type(OFile)  :: logfile
  type(String) :: temperature_dir
  type(OFile)  :: temperature_file
  type(OFile)  :: dispersion_file
  type(OFile)  :: symmetry_points_file
  type(OFile)  :: sampled_qpoints_file
  type(String) :: subspace_dir
  type(OFile)  :: vscha_castep_phonon_file
  type(OFile)  :: vscha_qe_force_constants_file
  type(OFile)  :: wavefunctions_file
  type(OFile)  :: vscha_thermodynamic_file
  type(OFile)  :: vscha1_thermodynamic_file
  type(OFile)  :: vscha2_thermodynamic_file
  type(OFile)  :: vscf_thermodynamic_file
  type(OFile)  :: pdos_file
  
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
  split_qpoints = lgcl(arguments%value('split_q-points'))
  energy_convergence = dble(arguments%value('energy_convergence'))
  no_converged_calculations_vscf = int(                  &
     & arguments%value('no_converged_calculations_vscf') )
  max_pulay_iterations = int(arguments%value('max_pulay_iterations'))
  pre_pulay_iterations = int(arguments%value('pre_pulay_iterations'))
  pre_pulay_damping = dble(arguments%value('pre_pulay_damping'))
  min_temperature = dble(arguments%value('min_temperature'))
  max_temperature = dble(arguments%value('max_temperature'))
  no_temperature_steps = int(arguments%value('no_temperature_steps'))
  min_frequency = dble(arguments%value('min_frequency'))
  path = arguments%value('path')
  no_path_points = int(arguments%value('no_path_points'))
  no_dos_samples = int(arguments%value('no_dos_samples'))
  no_vscf_basis_states = int(arguments%value('no_vscf_basis_states'))
  
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
  elseif (pre_pulay_damping<0_dp .or. pre_pulay_damping>1_dp) then
    call print_line(ERROR//': pre_pulay_damping must be between 0 and 1.')
    call quit()
  elseif (pre_pulay_iterations<2) then
    call print_line(ERROR//': pre_pulay_iterations must be at least 2.')
    call quit()
  elseif (max_pulay_iterations<2) then
    call print_line(ERROR//': max_pulay_iterations must be at least 2.')
    call quit()
  endif
  
  ! Read in setup_harmonic settings.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! Read in setup_anharmonic settings.
  setup_anharmonic_arguments = Dictionary(CaesarMode('setup_anharmonic'))
  call setup_anharmonic_arguments%read_file('setup_anharmonic.used_settings')
  calculate_stress = lgcl(setup_anharmonic_arguments%value('calculate_stress'))
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  modes = anharmonic_data%complex_modes
  qpoints = anharmonic_data%qpoints
  subspaces = anharmonic_data%degenerate_subspaces
  supercell = anharmonic_data%anharmonic_supercell
  
  ! Read in anharmonic potential.
  potential_file = IFile('potential.dat')
  potential = PotentialPointer(potential_file%lines())
  
  ! Read in anharmonic stress.
  if (calculate_stress) then
    stress_file = IFile('stress.dat')
    stress = StressPointer(stress_file%lines())
  endif
  
  ! --------------------------------------------------
  ! Generate objects which don't depend on the potential:
  !    - thermal energies from temperature range.
  !    - minimum image data.
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
  
  ! Make output directory, write out temperatures, and open logfile.
  output_dir = 'anharmonic_observables'
  call mkdir(output_dir)
  
  temperature_file = OFile(output_dir//'/temperatures.dat')
  call temperature_file%print_line('Thermal energies, KbT, (Ha):')
  call temperature_file%print_lines(thermal_energies)
  
  logfile = OFile(output_dir//'/anharmonic_observables_log.dat')
  
  ! Construct minimum image data.
  min_images = calculate_min_images(supercell)
  
  ! Generate stress prefactors, which are geometric factors quantifying
  !    the direction of each pair of modes w/r/t the lattice vectors.
  if (calculate_stress) then
    stress_prefactors = [( StressPrefactors(subspaces(i),modes), &
                         & i=1,                                  &
                         & size(subspaces)                       )]
  endif
  
  ! Calculate starting frequencies.
  starting_frequencies = [                                  &
     & max( subspaces%frequency,                            &
     &      anharmonic_data%frequency_of_max_displacement ) ]
  
  ! --------------------------------------------------
  ! Run calculations at each temperature.
  ! --------------------------------------------------
  allocate( vscha_frequencies( size(subspaces),                 &
          &                    size(thermal_energies) ),        &
          & vscf_basis(size(subspaces)),                        &
          & vscha_modes( anharmonic_data%structure%no_modes,    &
          &              size(qpoints)                       ), &
          & dynamical_matrices(size(qpoints)),                  &
          & vscf_thermodynamics( size(subspaces),               &
          &                      size(thermal_energies) ),      &
          & vscha_thermodynamics(size(thermal_energies)),       &
          & vscha1_thermodynamics( size(subspaces),             &
          &                        size(thermal_energies) ),    &
          & vscha2_thermodynamics( size(subspaces),             &
          &                        size(thermal_energies) ),    &
          & stat=ialloc); call err(ialloc)
  do i=1,size(thermal_energies)
    call print_line('')
    call print_line(repeat('=',50))
    call print_line('Thermal energy '//i//' of '//size(thermal_energies)// &
       & ': KT = '//thermal_energies(i)//' Ha')
    call print_line(repeat('=',50))
    
    ! Make temperature directory.
    temperature_dir = output_dir//'/temperature_'// &
       & left_pad(i,str(size(thermal_energies)))
    call mkdir(temperature_dir)
    
    ! --------------------------------------------------
    ! Run VSCF using VSCHA basis.
    ! --------------------------------------------------
    
    if (i==1) then
      vscha_basis = HarmonicBasis( subspaces%id,        &
                                 & starting_frequencies )
    else
      vscha_basis = HarmonicBasis( subspaces%id,            &
                                 & vscha_frequencies(:,i-1) )
    endif
    
    call print_line('Running VSCHA.')
    ! N.B. for the first run, vscha_output is not allocated,
    !    so starting_configuration will not be present() inside run_vscf.
    ! For later runs, starting_configuration will be the output of the
    !    previous run.
    vscha_output = run_vscf( potential,                          &
                           & stress,                             &
                           & subspaces,                          &
                           & vscha_basis,                        &
                           & thermal_energies(i),                &
                           & energy_convergence,                 &
                           & starting_frequencies,               &
                           & no_converged_calculations_vscf,     &
                           & max_pulay_iterations,               &
                           & pre_pulay_iterations,               &
                           & pre_pulay_damping,                  &
                           & anharmonic_data,                    &
                           & random_generator,                   &
                           & starting_configuration=vscha_output )
    vscha_potentials = [(vscha_output(j)%potential, j=1, size(vscha_output))]
    vscha_states = [( HarmonicStates(vscha_output(i)%states), &
                    & i=1,                                    &
                    & size(vscha_output)                      )]
    vscha_frequencies(:,i) = vscha_states%frequency
    if (calculate_stress) then
      vscha_stresses = [(vscha_output(j)%stress, j=1, size(vscha_output))]
    endif
    call print_line('')
    
    ! --------------------------------------------------
    ! Calculate thermodynamic quantities with VSCHA basis.
    ! --------------------------------------------------
    
    call print_line('Calculating thermodynamic observables under VSCHA.')
    do j=1,size(subspaces)
      if (calculate_stress) then
        subspace_stress = vscha_stresses(j)
        stress_prefactor = stress_prefactors(j)%average_prefactor()
      endif
      
      vscha1_thermodynamics(j,i) = harmonic_observables( &
         &    thermal_energy   = thermal_energies(i),    &
         &    stress           = subspace_stress,        &
         &    stress_prefactor = stress_prefactor,       &
         &    frequency        = vscha_frequencies(j,i), &
         &    num_dimensions   = size(subspaces(j)),     &
         &    anharmonic_data  = anharmonic_data)        &
         & / supercell%sc_size
      
      vscha2_thermodynamics(j,i) = effective_harmonic_observables( &
                   &    thermal_energy   = thermal_energies(i),    &
                   &    potential        = vscha_potentials(j),    &
                   &    stress           = subspace_stress,        &
                   &    stress_prefactor = stress_prefactor,       &
                   &    frequency        = vscha_frequencies(j,i), &
                   &    num_dimensions   = size(subspaces(j)),     &
                   &    anharmonic_data  = anharmonic_data)        &
                   & / supercell%sc_size
    enddo
    
    ! --------------------------------------------------
    ! Run VSCF using VCI basis.
    ! --------------------------------------------------
    
    ! Generate basis of states.
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
              & supercell,                                      &
              & no_vscf_basis_states-1,                         &
              & anharmonic_data%potential_expansion_order,      &
              & anharmonic_data%structure%symmetries            ))
      else
        vscf_basis(j) = SubspaceBasisPointer(FullSubspaceBasis( &
              & subspaces(j),                                   &
              & vscha_frequencies(j,i),                         &
              & modes,                                          &
              & qpoints,                                        &
              & supercell,                                      &
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
    vscf_output = run_vscf( potential,                         &
                          & stress,                            &
                          & subspaces,                         &
                          & vscf_basis,                        &
                          & thermal_energies(i),               &
                          & energy_convergence,                &
                          & vscha_frequencies(:,i),            &
                          & no_converged_calculations_vscf,    &
                          & max_pulay_iterations,              &
                          & pre_pulay_iterations,              &
                          & pre_pulay_damping,                 &
                          & anharmonic_data,                   &
                          & random_generator,                  &
                          & starting_configuration=vscf_output )
    
    subspace_potentials = [(vscf_output(j)%potential, j=1, size(vscf_output))]
    subspace_states = [(vscf_output(j)%states, j=1, size(vscf_output))]
    if (calculate_stress) then
      subspace_stresses = [(vscf_output(j)%stress, j=1, size(vscf_output))]
    endif
    call print_line('')
    
    ! --------------------------------------------------
    ! Calculate thermodynamic quantities with VCI basis.
    ! --------------------------------------------------
    
    call print_line('Calculating thermodynamic observables under VSCF.')
    if (calculate_stress) then
      vscf_thermodynamics(:,i) = vscf_basis%thermodynamic_data( &
                             &           thermal_energies(i),   &
                             &           subspace_states,       &
                             &           subspaces,             &
                             &           subspace_potentials,   &
                             &           subspace_stresses,     &
                             &           stress_prefactors,     &
                             &           anharmonic_data      ) &
                             & / supercell%sc_size
    else
      vscf_thermodynamics(:,i) = vscf_basis%thermodynamic_data(              &
                             &           thermal_energies(i),                &
                             &           subspace_states,                    &
                             &           subspaces,                          &
                             &           subspace_potentials,                &
                             &           anharmonic_data = anharmonic_data ) &
                             & / supercell%sc_size
    endif
    
    !! Print VSCF spectra and wavefunction information.
    !subspace_wavefunctions = SubspaceWavefunctionsPointer( &
    !         & vscf_basis%wavefunctions( subspace_states,  &
    !         &                           subspaces,        &
    !         &                           anharmonic_data ) )
    !do j=1,size(subspaces)
    !  subspace_dir = temperature_dir//'/subspace_'// &
    !     & left_pad(subspaces(j)%id, str(maxval(subspaces%id)))
    !  call mkdir(subspace_dir)
    !  
    !  wavefunctions_file = OFile(subspace_dir//'/vscf_wavefunctions.dat')
    !  call wavefunctions_file%print_lines(subspace_wavefunctions(j))
    !enddo
    
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
    hessian = reconstruct_hessian( supercell,          &
                                 & qpoints,            &
                                 & dynamical_matrices, &
                                 & logfile             )
    
    ! Write out Castep .phonon file and QE .fc file for the VSCHA modes.
    vscha_castep_phonon_file = OFile(                                &
       & temperature_dir//'/'//make_castep_phonon_filename(seedname) )
    call write_castep_phonon_file( vscha_castep_phonon_file, &
                                 & vscha_modes,              &
                                 & qpoints,                  &
                                 & anharmonic_data%structure )
    
    vscha_qe_force_constants_file = OFile(                                &
       & temperature_dir//'/'//make_qe_force_constants_filename(seedname) )
    call write_qe_force_constants_file( vscha_qe_force_constants_file, &
                                      & hessian,                       &
                                      & anharmonic_data%structure,     &
                                      & supercell                      )
    
    ! Generate self-consistent phonon dispersion curve by interpolating between
    !    calculated q-points using Fourier interpolation.
    phonon_dispersion = PhononDispersion( supercell,      &
                                        & min_images,     &
                                        & hessian,        &
                                        & path,           &
                                        & no_path_points, &
                                        & logfile         )
    
    symmetry_points_file = OFile(temperature_dir//'/high_symmetry_points.dat')
    call symmetry_points_file%print_lines(phonon_dispersion%path())
    
    dispersion_file = OFile(temperature_dir//'/phonon_dispersion_curve.dat')
    call dispersion_file%print_lines(phonon_dispersion%frequencies())
    
    ! Generate self-consistent phonon density of states,
    !    interpolating as above.
    ! Re-seed the random generator each time, so that the results at different
    !    temperatures can be compared free from random noise.
    random_generator = RandomReal(random_generator_seed)
    phonon_dos = PhononDos( supercell,             &
                          & min_images,            &
                          & hessian,               &
                          & [thermal_energies(i)], &
                          & min_frequency,         &
                          & no_dos_samples,        &
                          & logfile,               &
                          & random_generator       )
    
    ! Write out dos and dispersion at this temperature.
    sampled_qpoints_file = OFile(temperature_dir//'/sampled_qpoints.dat')
    call sampled_qpoints_file%print_line('q-point (x,y,z) | &
                                         &number of frequencies ignored')
    call sampled_qpoints_file%print_lines(phonon_dos%qpoints)
    
    pdos_file = OFile(temperature_dir//'/phonon_density_of_states.dat')
    call pdos_file%print_lines(phonon_dos%pdos)
    
    vscha_thermodynamics(i) = phonon_dos%thermodynamic_data(1)
    
    ! Interpolate stress under the effective harmonic approximation.
    if (calculate_stress) then
      if (stress%can_be_interpolated()) then
        call vscha_thermodynamics(i)%set_stress(      &
           & stress%calculate_interpolated_stress(    &
           &            1e-30_dp,                     &
           &            phonon_dos%qpoints%qpoint,    &
           &            thermal_energies(i),          &
           &            min_frequency,                &
           &            supercell,                    &
           &            hessian,                      &
           &            min_images,                   &
           &            subspaces,                    &
           &            vscha_basis,                  &
           &            vscha_states,                 &
           &            min_images,                   &
           &            anharmonic_data            ), &
           & anharmonic_data%structure%volume         )
      endif
    endif
    
    ! --------------------------------------------------
    ! Calculated interpolated VSCF results.
    ! --------------------------------------------------
    if (potential%can_be_interpolated()) then
      call print_line('Interpolating VSCF results.')
      call potential%calculate_interpolated_thermodynamics( &
                                     & thermal_energies(i), &
                                     & min_frequency,       &
                                     & subspaces,           &
                                     & subspace_potentials, &
                                     & vscf_basis,          &
                                     & subspace_states,     &
                                     & anharmonic_data      )
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
  vscha1_thermodynamic_file = OFile(                    &
     & output_dir//'/vscha_thermodynamic_variables.dat' )
  call vscha1_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscha1_thermodynamic_file%print_line(sum(vscha1_thermodynamics(:,i)))
  enddo
  
  ! VSCHA without Fourier interpolation, but with <V> taken using the VSCF
  !    potential rather than the effective harmonic potential.
  vscha2_thermodynamic_file = OFile(                         &
     & output_dir//'/vscha_vscf_thermodynamic_variables.dat' )
  call vscha2_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscha2_thermodynamic_file%print_line(sum(vscha2_thermodynamics(:,i)))
  enddo
  
  ! VSCF without Fourier interpolation.
  vscf_thermodynamic_file = OFile(                     &
     & output_dir//'/vscf_thermodynamic_variables.dat' )
  call vscf_thermodynamic_file%print_line(thermodynamic_header)
  do i=1,size(thermal_energies)
    call vscf_thermodynamic_file%print_line(sum(vscf_thermodynamics(:,i)))
  enddo
  
  ! VSCHA with Fourier interpolation.
  vscha_thermodynamic_file = OFile(                                  &
     & output_dir//'/interpolated_vscha_thermodynamic_variables.dat' )
  call vscha_thermodynamic_file%print_line(thermodynamic_header)
  call vscha_thermodynamic_file%print_lines(vscha_thermodynamics)
end subroutine
end module
