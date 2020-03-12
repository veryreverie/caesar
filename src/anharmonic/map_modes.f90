! ======================================================================
! Maps the potential along each normal mode.
! ======================================================================
module map_modes_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  
  use mode_map_module
  implicit none
  
  private
  
  public :: startup_map_modes
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_map_modes()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'map_modes'
  mode%description = 'Maps the potential along normal modes. If use_potential &
  &is true, map_modes should be run after calculate_potential, &
  &otherwise it may be run after calculate_normal_modes.'
  mode%keywords = [                                                           &
  & KeywordData( 'harmonic_path',                                             &
  &              'harmonic_path is the path to the directory where harmonic &
  &calculations were run.',                                                   &
  &              default_value='.',                                           &
  &              is_path=.true.),                                             &
  & KeywordData( 'q-point_grid',                                              &
  &              'q-point_grid is the number of q-points in each direction in &
  &a Monkhorst-Pack grid. This should be specified as three integers &
  &separated by spaces. All q-points in this grid must also appear in the &
  &grid used for harmonic calculations.'),                                    &
  & KeywordData( 'mode_ids',                                                  &
  &              'mode_ids is the list of mode IDs of the real modes which &
  &should be mapped. If mode_ids is not present then all modes will be &
  &mapped. mode_ids should be given as space-separated integers.',            &
  &              is_optional=.true.),                                         &
  & KeywordData( 'no_single_mode_samples',                                    &
  &              'no_single_mode_samples is the number of points (either &
  &side of zero) along each mode at which the anharmonic potential will be &
  &sampled.',                                                                 &
  &              default_value='100'),                                        &
  & KeywordData( 'use_potential',                                             &
  &              'use_potential specifies whether or not to use the &
  &anharmonic potential. The anharmonic potential must have been &
  &calculated using calculate_potential in order for this to be true.',       &
  &              default_value='false'),                                      &
  & KeywordData( 'maximum_displacement',                                      &
  &              'maximum_displacement, u_max, is the largest distance any &
  &sampling point will be from the equilibrium position. maximum_displacement &
  &should be given in Bohr. Due to the use of mass-reduced co-ordinates, in &
  &systems containing different elements modes with higher contributions from &
  &heavier atoms will be displaced less far than this, by a factor of &
  &sqrt(m_min/m), where m is the mass of the element and m_min is the minimum &
  &atomic mass in the structure.'),                                           &
  & KeywordData( 'max_energy_of_displacement',                                &
  &              'max_energy_of_displacement, E_max is the energy of the &
  &harmonic potential up to which the anharmonic potential is sampled. A mode &
  &with harmonic frequency w and effective mass m will be displaced up to &
  &0.5*m*(wu)^2=E_max, unless this would lead to u>u_max. This is equivalent &
  &to E_max = 0.5 * m_min * (w_max*u_max)^2. E_max should be given in &
  &Hartree.',                                                                 &
  &              exclusive_with=[str('frequency_of_max_displacement')]),      &
  & KeywordData( 'frequency_of_max_displacement',                             &
  &              'frequency_of_max_displacement is the frequency, w_min, at &
  &which maximum displacement happens. Displacement along modes with w>w_min &
  &is scaled by sqrt(w_min/w), and displacement along modes with w<w_min &
  &is unscaled. This is equivalent to E_max = 0.5 * m_min * (w_max*u_max)^2. &
  &w_min should be given in Hartree.',                                        &
  &              exclusive_with=[str('max_energy_of_displacement')]),         &
  & KeywordData( 'calculate_stress',                                          &
  &              'calculate_stress specifies whether or not to calculate &
  &stress and pressure. If calculate_stress is true then calculate_stress &
  &for setup_anharmonic must also have been true.',                           &
  &              default_value='true'),                                       &
  & KeywordData( 'validate_potential',                                        &
  &              'validate_potential specifies that the potential should &
  &be verified against fresh electronic structure calculations when &
  &sampling. Depending on the electronic structure method, this is likely &
  &to be very computationally intensive.',                                    &
  &              default_value='false'),                                      &
  & KeywordData( 'run_script',                                                &
  &              'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.),                                            &
  & KeywordData( 'no_cores',                                                  &
  &              'no_cores is the number of cores on which the electronic &
  &structure calculation will be run. This is passed to the specified run &
  &script.',                                                                  &
  &              default_value='1'),                                          &
  & KeywordData( 'no_nodes',                                                  &
  &              'no_nodes is the number of nodes on which the electronic &
  &structure calculation will be run. This is passed to the specified run &
  &script.',                                                                  &
  &              default_value='1'),                                          &
  & KeywordData( 'run_script_data',                                           &
  &              'run_script_data will be passed to the specified run &
  &script after all other arguments. This should be used to pass &
  &information not covered by the other arguments.',                          &
  &              default_value=''),                                           &
  & KeywordData( 'calculation_type',                                          &
  &              'calculation_type specifies whether any electronic &
  &structure calculations should be run in addition to the user-defined &
  &script. Settings are: "none" and "quip".',                                 &
  &              default_value='none')                                        ]
  mode%main_subroutine => map_modes_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine map_modes_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  type(String)          :: harmonic_path
  integer               :: qpoint_grid(3)
  integer,  allocatable :: mode_ids(:)
  integer               :: no_single_mode_samples
  logical               :: use_potential
  logical               :: calculate_stress
  logical               :: validate_potential
  real(dp)              :: maximum_displacement
  real(dp), allocatable :: max_energy_of_displacement
  real(dp), allocatable :: frequency_of_max_displacement
  type(String)          :: run_script
  integer               :: no_cores
  integer               :: no_nodes
  type(String)          :: run_script_data
  type(String)          :: calculation_type
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Harmonic data.
  type(StructureData)                :: structure
  type(QpointData),      allocatable :: harmonic_qpoints(:)
  type(DynamicalMatrix), allocatable :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     allocatable :: harmonic_complex_modes(:,:)
  
  ! Max displacement data.
  type(MaxDisplacement) :: max_displacement
  
  ! Anharmonic data.
  type(AnharmonicData)                  :: anharmonic_data
  type(StructureData)                   :: large_supercell
  type(QpointData),         allocatable :: qpoints(:)
  type(RealMode),           allocatable :: real_modes(:)
  type(DegenerateSubspace), allocatable :: subspaces(:)
  
  ! Interpolated supercell.
  type(InterpolatedSupercell) :: interpolated_supercell
  
  ! Map from mode ids to locations in mode_ids.
  integer, allocatable :: mode_locs(:)
  integer              :: loc
  
  ! Anharmonic potential and stress.
  type(PotentialPointer), allocatable :: potential
  class(StressData),      allocatable :: stress
  
  ! Electronic structure calculation handlers.
  type(CalculationWriter) :: calculation_writer
  type(CalculationRunner) :: calculation_runner
  type(CalculationReader) :: calculation_reader
  
  ! Variables for sampline potential.
  real(dp),      allocatable :: unscaled_mode_displacements(:)
  real(dp),      allocatable :: mode_displacements(:)
  real(dp),      allocatable :: l2_cartesian_displacements(:)
  type(ModeMap), allocatable :: mode_maps(:)
  integer,       allocatable :: qpoint_modes(:)
  integer,       allocatable :: subspace_modes(:)
  
  ! Variables for validating potential.
  real(dp), allocatable       :: sampled_energies(:)
  type(RealModeForce)         :: sampled_force
  real(dp), allocatable       :: sampled_forces(:)
  real(dp), allocatable       :: sampled_pressures(:)
  type(IntMatrix)             :: supercell_matrix
  type(StructureData)         :: supercell
  type(RealMode)              :: mode
  type(RealModeDisplacement)  :: real_mode_displacement
  type(CartesianDisplacement) :: cartesian_displacement
  type(StructureData)         :: displaced_structure
  type(ElectronicStructure)   :: electronic_structure
  
  ! Files and directories.
  type(IFile)  :: structure_file
  type(IFile)  :: harmonic_qpoints_file
  type(IFile)  :: harmonic_dynamical_matrices_file
  type(IFile)  :: harmonic_complex_modes_file
  type(IFile)  :: anharmonic_data_file
  type(IFile)  :: potential_file
  type(IFile)  :: stress_file
  type(String) :: qpoint_dir
  type(String) :: subspace_dir
  type(String) :: output_dir
  type(OFile)  :: supercell_file
  type(OFile)  :: mode_maps_file
  type(String) :: mode_dir
  type(String) :: displacement_dir
  
  ! Temporary variables.
  integer :: i,j,k,ialloc
  
  ! --------------------------------------------------
  ! Read in inputs and previously calculated data.
  ! --------------------------------------------------
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split_line(arguments%value('q-point_grid')))
  if (arguments%is_set('mode_ids')) then
    mode_ids = int(tokens(arguments%value('mode_ids')))
  endif
  no_single_mode_samples = int(arguments%value('no_single_mode_samples'))
  use_potential = lgcl(arguments%value('use_potential'))
  calculate_stress = lgcl(arguments%value('calculate_stress'))
  validate_potential = lgcl(arguments%value('validate_potential'))
  maximum_displacement = dble(arguments%value('maximum_displacement'))
  if (arguments%is_set('frequency_of_max_displacement')) then
    frequency_of_max_displacement = &
       & dble(arguments%value('frequency_of_max_displacement'))
  else
    max_energy_of_displacement = &
       & dble(arguments%value('max_energy_of_displacement'))
  endif
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  no_nodes = int(arguments%value('no_nodes'))
  run_script_data = arguments%value('run_script_data')
  calculation_type = arguments%value('calculation_type')
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file(                 &
          & harmonic_path//'/setup_harmonic.used_settings' )
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  
  ! Read in harmonic data.
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  ! Construct max displacement data.
  max_displacement = MaxDisplacement(                                 &
     & maximum_displacement          = maximum_displacement,          &
     & structure                     = structure,                     &
     & frequency_of_max_displacement = frequency_of_max_displacement, &
     & max_energy_of_displacement    = max_energy_of_displacement     )
  
  if (use_potential) then
    ! Read in anharmonic data.
    anharmonic_data_file = IFile('anharmonic_data.dat')
    anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
    
    large_supercell = anharmonic_data%anharmonic_supercell
    qpoints = anharmonic_data%qpoints
    real_modes = anharmonic_data%real_modes
    subspaces = anharmonic_data%degenerate_subspaces
  else
    harmonic_qpoints_file = IFile(harmonic_path//'/qpoints.dat')
    harmonic_qpoints = QpointData(harmonic_qpoints_file%sections())
    
    allocate( harmonic_dynamical_matrices(size(harmonic_qpoints)), &
            & harmonic_complex_modes( structure%no_modes,          &
            &                         size(harmonic_qpoints) ),    &
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
    
    interpolated_supercell = InterpolatedSupercell( &
                     & qpoint_grid,                 &
                     & structure,                   &
                     & harmonic_qpoints,            &
                     & harmonic_dynamical_matrices, &
                     & harmonic_complex_modes       )
    
    large_supercell = interpolated_supercell%supercell
    qpoints         = interpolated_supercell%qpoints
    real_modes      = interpolated_supercell%real_modes
  endif
  
  ! Construct mode locations array.
  if (.not. arguments%is_set('mode_ids')) then
    mode_ids = real_modes%id
  endif
  mode_locs = [(0,i=1,maxval(real_modes%id))]
  j = 0
  do i=1,size(real_modes)
    if (real_modes(i)%id .in. mode_ids) then
      j = j+1
      mode_locs(real_modes(i)%id) = j
    endif
  enddo
  
  ! Read in anharmonic potential and stress.
  if (use_potential) then
    potential_file = IFile('potential.dat')
    potential = PotentialPointer(potential_file%lines())
    if (calculate_stress) then
      stress_file = IFile('stress.dat')
      stress = StressPointer(stress_file%lines())
    endif
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation handlers.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  calculation_runner = CalculationRunner(      &
     & file_type           = file_type,        &
     & seedname            = seedname,         &
     & run_script          = run_script,       &
     & no_cores            = no_cores,         &
     & no_nodes            = no_nodes,         &
     & run_script_data     = run_script_data,  &
     & calculation_type    = calculation_type, &
     & use_forces          = .true.,           &
     & use_hessians        = .false.,          &
     & calculate_stress    = calculate_stress, &
     & exit_on_error       = .true.,           &
     & repeat_calculations = .true.            )
  
  calculation_reader = CalculationReader()
  
  ! --------------------------------------------------
  ! Map the potential.
  ! --------------------------------------------------
  ! Calculate displacements before scaling by 1/sqrt(frequency).
  unscaled_mode_displacements =                                 &
     &   [(i,i=-no_single_mode_samples,no_single_mode_samples)] &
     & * max_displacement%maximum_weighted_displacement         &
     & / no_single_mode_samples
  
  allocate( mode_maps(size(real_modes)),                                   &
          & l2_cartesian_displacements(size(unscaled_mode_displacements)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(real_modes)
    loc = mode_locs(real_modes(i)%id)
    if (loc/=0) then
      ! Scale displacement by 1/sqrt(frequency).
      mode_displacements =                                               &
         &   unscaled_mode_displacements                                 &
         & * sqrt( max_displacement%frequency_of_max_displacement        &
         &       / max( real_modes(i)%frequency,                         &
         &              max_displacement%frequency_of_max_displacement ) )
      
      ! Calculate the L2 displacments in Bohr.
      do j=1,size(mode_displacements)
        real_mode_displacement = RealModeDisplacement( &
                             & [real_modes(i)],        &
                             & [mode_displacements(j)] )
        cartesian_displacement = CartesianDisplacement( &
                              & real_mode_displacement, &
                              & large_supercell,        &
                              & real_modes,             &
                              & qpoints                 )
        l2_cartesian_displacements(j) =                  &
           & sqrt(sum( cartesian_displacement%vectors    &
           &         * cartesian_displacement%vectors )) &
           & / large_supercell%sc_size
      enddo
      l2_cartesian_displacements(:no_single_mode_samples) = &
         & -l2_cartesian_displacements(:no_single_mode_samples)
      
      ! Sample the model potential.
      if (use_potential) then
        mode_maps(loc) = ModeMap( mode_displacements,         &
                                & l2_cartesian_displacements, &
                                & real_modes(i),              &
                                & potential,                  &
                                & stress,                     &
                                & anharmonic_data             )
      else
        mode_maps(loc) = ModeMap( mode_displacements,         &
                                & l2_cartesian_displacements, &
                                & real_modes(i)               )
      endif
    endif
  enddo
  
  ! --------------------------------------------------
  ! Validate the model potential by sampling the electronic structure at
  !    the same points as the model potential.
  ! --------------------------------------------------
  allocate( sampled_energies(size(mode_displacements)), &
          & sampled_forces(size(mode_displacements)),   &
          & stat=ialloc); call err(ialloc)
  if (calculate_stress) then
    allocate( sampled_pressures(size(mode_displacements)), &
            & stat=ialloc); call err(ialloc)
  endif
  do i=1,size(qpoints)
    if (qpoints(i)%paired_qpoint_id<qpoints(i)%id) then
      cycle
    endif
    qpoint_modes = filter(real_modes%qpoint_id_plus==qpoints(i)%id)
    if (qpoints(i)%paired_qpoint_id>qpoints(i)%id) then
      qpoint_modes = [ qpoint_modes,                                     &
                     & filter(real_modes%qpoint_id_minus==qpoints(i)%id) ]
    endif
    
    if (all(mode_locs(real_modes(qpoint_modes)%id)==0)) then
      cycle
    endif
    
    qpoint_dir = 'qpoint_'//left_pad(qpoints(i)%id, str(maxval(qpoints%id)))
    call mkdir(qpoint_dir)
    
    ! Construct and write out supercell.
    supercell_matrix = construct_supercell_matrix(qpoints(i), structure)
    supercell = construct_supercell( structure,       &
                                   & supercell_matrix )
    supercell_file = OFile(qpoint_dir//'/structure.dat')
    call supercell_file%print_lines(supercell)
    
    do j=1,size(qpoint_modes)
      mode = real_modes(qpoint_modes(j))
      loc = mode_locs(mode%id)
      if (loc/=0) then
        mode_dir = qpoint_dir//'/real_mode_'// &
                 & left_pad(mode%id, str(maxval(real_modes%id)))
        call mkdir(mode_dir)
        do k=1,size(sampled_energies)
          real_mode_displacement = RealModeDisplacement( &
                & [mode],                                &
                & [mode_maps(loc)%mode_displacements(k)] )
          cartesian_displacement = CartesianDisplacement( &
                                & real_mode_displacement, &
                                & supercell,              &
                                & real_modes,             &
                                & qpoints                 )
          displaced_structure = displace_structure( supercell,             &
                                                  & cartesian_displacement )
          
          displacement_dir = mode_dir//'/displacement_'// &
                           & left_pad(k,str(size(mode_displacements)))
          call calculation_writer%write_calculation( &
                              & displaced_structure, &
                              & displacement_dir     )
          
          if (validate_potential) then
            call calculation_runner%run_calculation(displacement_dir)
            
            electronic_structure = calculation_reader%read_calculation( &
                                                     & displacement_dir )
            
            sampled_energies(k) = electronic_structure%energy()
            sampled_force = RealModeForce( electronic_structure%forces(), &
                                         & supercell,                     &
                                         & real_modes,                    &
                                         & qpoints                        )
            sampled_forces(k) = sampled_force%force(mode)
            if (calculate_stress) then
              sampled_pressures(k) = trace(electronic_structure%stress())/3
            endif
          endif
        enddo
        
        if (validate_potential) then
          sampled_energies = ( sampled_energies                             &
                         &   - sampled_energies(no_single_mode_samples+1) ) &
                         & / supercell%sc_size
          if (calculate_stress) then
            sampled_pressures = sampled_pressures &
                            & - sampled_pressures(no_single_mode_samples+1)
          endif
          mode_maps(loc)%sampled_energies = sampled_energies
          mode_maps(loc)%sampled_forces   = sampled_forces
          if (calculate_stress) then
            mode_maps(loc)%sampled_pressures = sampled_pressures
          endif
        endif
      endif
    enddo
  enddo
  
  ! --------------------------------------------------
  ! Write out output.
  ! --------------------------------------------------
  if (use_potential) then
    do i=1,size(subspaces)
      subspace_modes = filter(real_modes%id .in. subspaces(i)%mode_ids)
      if (any(mode_locs(real_modes(subspace_modes)%id)/=0)) then
        subspace_dir = 'subspace_'//left_pad( subspaces(i)%id,          &
                                            & str(maxval(subspaces%id)) )
        call mkdir(subspace_dir)
        mode_maps_file = OFile(subspace_dir//'/mode_maps.dat')
        call mode_maps_file%print_line(                    &
           & 'Harmonic frequencies: '//subspaces%frequency )
        do j=1,size(subspace_modes)
          loc = mode_locs(real_modes(subspace_modes(j))%id)
          if (loc/=0) then
            call mode_maps_file%print_line('')
            call mode_maps_file%print_lines(mode_maps(loc))
          endif
        enddo
      endif
    enddo
  else
    output_dir = 'mode_maps'
    call mkdir(output_dir)
    mode_maps_file = OFile(output_dir//'/mode_maps.dat')
    do i=1,size(real_modes)
      loc = mode_locs(real_modes(i)%id)
      if (loc/=0) then
        call mode_maps_file%print_line('')
        call mode_maps_file%print_lines(mode_maps(loc))
      endif
    enddo
  endif
end subroutine
end module
