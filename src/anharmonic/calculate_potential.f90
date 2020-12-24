! ======================================================================
! Calculates the anharmonic potential, using the results of run_anharmonic.
! ======================================================================
module calculate_potential_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  use interpolation_module
  implicit none
  
  private
  
  public :: startup_calculate_potential
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_calculate_potential()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'calculate_potential'
  mode%description = 'Uses the results of run_anharmonic to calculate &
     &the anharmonic potential. Should be run after run_anharmonic.'
  mode%keywords = [                                                           &
     & KeywordData( 'energy_to_force_ratio',                                  &
     &              'energy_to_force_ratio is the same as &
     &energy_to_force_ratio in setup_anharmonic. If unset, this will default &
     &to the value set in setup_anharmonic. This may be varied here to check &
     &that it does not have a large effect on the observables.',              &
     &              is_optional=.true.),                                      &
     & KeywordData( 'interpolate_potential',                                  &
     &              'interpolate_potential determines whether or not the &
     &potential should be interpolated before vscf and statistical mechanical &
     &calculations are performed. The interpolated_ keywords are only used if &
     &interpolate_potential is true.',                                        &
     &              default_value='false'),                                   &
     & KeywordData( 'interpolated_q-point_grid',                              &
     &              'interpolated_q-point_grid is the number of q-points in &
     &each direction in the Monkhorst-Pack grid onto which the potential will &
     &be interpolated. This should be specified as three integers separated &
     &by spaces.'),                                                           &
     & KeywordData( 'interpolated_maximum_coupling_order',                    &
     &              'interpolated_maximum_coupling_order is the maximum &
     &number of degenerate subspaces which may be coupled together in the &
     &interpolated potential. Must be at least 1.'),                          &
     & KeywordData( 'interpolated_potential_expansion_order',                 &
     &              'interpolated_potential_expansion_order is the order up &
     &to which the interpolated potential is expanded. e.g. if &
     &interpolated_potential_expansion_order=4 then terms up to and including &
     &u^4 are included. Must be at least 2, and at least as large as &
     &interpolated_maximum_coupling_order, and no larger than &
     &potential_expansion_order from setup_anharmonic.'),                     &
     & KeywordData( 'interpolated_vscf_basis_functions_only',                 &
     &              'interpolated_vscf_basis_functions_only specifies that &
     &the interpolated potential will only be expanded in terms of basis &
     &functions which are relevant to vscf.',                                 &
     &              default_value='true'),                                    &
     & KeywordData( 'loto_direction',                                         &
     &              'loto_direction specifies the direction (in reciprocal &
     &co-ordinates from which the gamma point is approached when calculating &
     &LO/TO corrections. See setup_harmonic help for details. loto_direction &
     &may only be specified here if it was not specified in setup_harmonic, &
     &and if every symmetry of the system leaves it invariant, i.e. q.S=q for &
     &all S, where S is the symmetry matrix and q is loto_direction. See &
     &structure.dat for the list of symmetries.',                             &
     &              is_optional = .true.)                                     ]
  mode%main_subroutine => calculate_potential_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine calculate_potential_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Input arguments.
  real(dp)             :: weighted_energy_force_ratio
  logical              :: interpolate_potential
  integer              :: interpolated_qpoint_grid(3)
  integer              :: interpolated_coupling_order
  integer              :: interpolated_expansion_order
  logical              :: interpolated_vscf_only
  type(Fractionvector) :: loto_direction
  logical              :: loto_direction_set
  
  ! Arguments to setup_harmonic.
  type(Dictionary) :: setup_harmonic_arguments
  
  ! Arguments to setup_anharmonic.
  type(Dictionary) :: setup_anharmonic_arguments
  type(String)     :: potential_representation
  integer          :: potential_expansion_order
  real(dp)         :: energy_to_force_ratio
  logical          :: calculate_stress
  
  ! Electronic structure calculation reader.
  type(CalculationReader) :: calculation_reader
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! Anharmonic potential.
  type(PotentialPointer) :: potential
  
  ! Harmonic variables.
  type(StructureData)                :: harmonic_supercell
  type(QpointData),      allocatable :: harmonic_qpoints(:)
  type(DynamicalMatrix), allocatable :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     allocatable :: harmonic_complex_modes(:,:)
  type(MinImages),       allocatable :: harmonic_min_images(:,:)
  type(DynamicalMatrix), allocatable :: difference_dynamical_matrices(:)
  
  ! Variables for interpolating the potential.
  type(MinImages), allocatable :: anharmonic_min_images(:,:)
  type(InterpolatedSupercell)  :: interpolated_supercell
  type(MaxDisplacement)        :: max_displacement
  type(AnharmonicData)         :: interpolated_anharmonic_data
  type(PotentialPointer)       :: interpolated_potential
  
  ! Stress
  type(SubspaceCoupling), allocatable :: stress_subspace_coupling(:)
  type(StressPointer)                 :: stress
  
  ! Files and directories.
  type(IFile)  :: anharmonic_data_file
  type(String) :: sampling_points_dir
  type(OFile)  :: logfile
  type(OFile)  :: potential_file
  type(String) :: harmonic_path
  type(IFile)  :: harmonic_supercell_file
  type(IFile)  :: harmonic_qpoints_file
  type(String) :: qpoint_dir
  type(IFile)  :: harmonic_dynamical_matrices_file
  type(IFile)  :: harmonic_complex_modes_file
  type(OFile)  :: interpolated_anharmonic_data_file
  type(OFile)  :: interpolated_potential_file
  type(OFile)  :: stress_subspace_coupling_file
  type(OFile)  :: stress_file
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! Parse input arguments.
  interpolate_potential = lgcl(arguments%value('interpolate_potential'))
  if (interpolate_potential) then
    interpolated_qpoint_grid = int(split_line(        &
       & arguments%value('interpolated_q-point_grid') ))
    interpolated_coupling_order = int(                          &
       & arguments%value('interpolated_maximum_coupling_order') )
    interpolated_expansion_order = int(                            &
       & arguments%value('interpolated_potential_expansion_order') )
    interpolated_vscf_only = lgcl(                                 &
       & arguments%value('interpolated_vscf_basis_functions_only') )
  endif
  
  ! Read in setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  
  ! Read in setup_anharmonic arguments.
  setup_anharmonic_arguments = Dictionary(CaesarMode('setup_anharmonic'))
  call setup_anharmonic_arguments%read_file( &
          & 'setup_anharmonic.used_settings' )
  potential_representation = &
     & setup_anharmonic_arguments%value('potential_representation')
  potential_expansion_order = &
     & int(setup_anharmonic_arguments%value('potential_expansion_order'))
  calculate_stress = lgcl(setup_anharmonic_arguments%value('calculate_stress'))
  harmonic_path = setup_anharmonic_arguments%value('harmonic_path')
  if (arguments%is_set('energy_to_force_ratio')) then
    energy_to_force_ratio = &
       & dble(arguments%value('energy_to_force_ratio'))
  else
    energy_to_force_ratio = &
       & dble(setup_anharmonic_arguments%value('energy_to_force_ratio'))
  endif
  
  ! Read in anharmonic data.
  anharmonic_data_file = IFile('anharmonic_data.dat')
  anharmonic_data = AnharmonicData(anharmonic_data_file%lines())
  
  ! Initialise LO/TO splitting if necessary.
  loto_direction_set = .false.
  if (setup_harmonic_arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
  endif
  
  if (arguments%is_set('loto_direction')) then
    if (loto_direction_set) then
      call print_line(ERROR//': loto_direction may not be specified here, &
         &since it was already specified in setup_harmonic.')
      call quit()
    endif
    loto_direction = FractionVector(                      &
       & setup_harmonic_arguments%value('loto_direction') )
    loto_direction_set = .true.
    if (any(loto_breaks_symmetry(                     &
       & anharmonic_data%structure%symmetries%tensor, &
       & loto_direction                               ))) then
      call print_line(ERROR//': loto_direction has been specified in a &
         &direction which breaks symmetry. To specify this direction, please &
         &set loto_direction when running setup_harmonic.')
      call quit()
    endif
  endif
  
  ! Initialise calculation reader for Potential mapping.
  if (loto_direction_set) then
    calculation_reader = CalculationReader(loto_direction)
  else
    calculation_reader = CalculationReader()
  endif
  
  ! Calculate weighted energy to force ratio.
  weighted_energy_force_ratio = &
     &   energy_to_force_ratio  &
     & * sqrt(maxval(anharmonic_data%structure%atoms%mass()))
  
  ! Initialise potential to the chosen representation
  if (potential_representation=='polynomial') then
    potential = PotentialPointer(PolynomialPotential(anharmonic_data))
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
  ! Open logfile.
  logfile = OFile('setup_anharmonic_logfile.dat')
  
  ! Generate the potential itself, and write it to file.
  call print_line('Generating potential energy surface.')
  sampling_points_dir = 'sampling_points'
  call potential%generate_potential( anharmonic_data,             &
                                   & weighted_energy_force_ratio, &
                                   & sampling_points_dir,         &
                                   & calculation_reader,          &
                                   & logfile                      )
  
  potential_file = OFile('potential.dat')
  call potential_file%print_lines(potential)
  
  ! Interpolate the potential.
  if (interpolate_potential) then
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
    
    harmonic_min_images = calculate_min_images(harmonic_supercell)
    
    anharmonic_min_images = calculate_min_images( &
           & anharmonic_data%anharmonic_supercell )
    
    interpolated_supercell = InterpolatedSupercell( &
                     & interpolated_qpoint_grid,    &
                     & anharmonic_data%structure,   &
                     & harmonic_supercell,          &
                     & harmonic_qpoints,            &
                     & harmonic_dynamical_matrices, &
                     & harmonic_complex_modes,      &
                     & logfile                      )
    
    ! max_displacement is not used, so a dummy is created.
    max_displacement = MaxDisplacement(                             &
       & maximum_displacement          = 0.0_dp,                    &
       & structure                     = anharmonic_data%structure, &
       & frequency_of_max_displacement = 0.0_dp,                    &
       & max_energy_of_displacement    = 0.0_dp                     )
    
    interpolated_anharmonic_data = AnharmonicData(                     &
       & structure                     = anharmonic_data%structure,    &
       & interpolated_supercell        = interpolated_supercell,       &
       & max_displacement              = max_displacement,             &
       & potential_expansion_order     = interpolated_expansion_order, &
       & maximum_coupling_order        = interpolated_coupling_order,  &
       & vscf_basis_functions_only     = interpolated_vscf_only,       &
       & energy_to_force_ratio         = 0.0_dp                        )
    
    interpolated_anharmonic_data_file = OFile( &
          & 'interpolated_anharmonic_data.dat' )
    call interpolated_anharmonic_data_file%print_lines( &
                         & interpolated_anharmonic_data )
    
    if (potential_representation=='polynomial') then
      interpolated_potential = PotentialPointer( &
          & PolynomialPotential(anharmonic_data) )
    else
      call print_line( ERROR//': Unrecognised potential representation: '// &
                     & potential_representation)
      call err()
    endif
    
    ! Calculate the harmonic terms which have longer range than the anharmonic
    !    supercell can represent.
    difference_dynamical_matrices = calculate_difference_dynamical_matrices( &
                                     & harmonic_dynamical_matrices,          &
                                     & harmonic_qpoints,                     &
                                     & harmonic_supercell,                   &
                                     & harmonic_min_images,                  &
                                     & anharmonic_data%qpoints,              &
                                     & anharmonic_data%anharmonic_supercell, &
                                     & anharmonic_min_images                 )
    difference_dynamical_matrices = interpolate_dynamical_matrices( &
                             & difference_dynamical_matrices,       &
                             & harmonic_qpoints,                    &
                             & harmonic_supercell,                  &
                             & harmonic_min_images,                 &
                             & interpolated_anharmonic_data%qpoints )
    
    ! Interpolate the potential, including the calculated anharmonic potential
    !    and the parts of the calculated harmonic potential which have longer
    !    range than the anharmonic supercell can represent.
    call interpolated_potential%interpolate_potential(                  &
       & anharmonic_min_images         = anharmonic_min_images,         &
       & potential                     = potential,                     &
       & anharmonic_data               = anharmonic_data,               &
       & interpolated_anharmonic_data  = interpolated_anharmonic_data,  &
       & difference_dynamical_matrices = difference_dynamical_matrices, &
       & logfile                       = logfile                        )
    
    interpolated_potential_file = OFile('interpolated_potential.dat')
    call interpolated_potential_file%print_lines(interpolated_potential)
  endif
  
  ! If calculate_stress is true, generate the stress and write it to file.
  if (calculate_stress) then
    call print_line('Generating stress surface.')
    ! Re-initialise calculation reader for stress mapping.
    if (loto_direction_set) then
      calculation_reader = CalculationReader(loto_direction)
    else
      calculation_reader = CalculationReader()
    endif
    
    ! Generates coupled subspaces for stress calculation.
    stress_subspace_coupling = generate_coupled_subspaces( &
                   & anharmonic_data%degenerate_subspaces, &
                   & 1                                     )
    stress_subspace_coupling_file = OFile('stress_subspace_coupling.dat')
    call stress_subspace_coupling_file%print_lines(stress_subspace_coupling)
    
    ! Generate stress.
    stress = potential%generate_stress(                        &
       & anharmonic_data           = anharmonic_data,          &
       & sampling_points_dir       = sampling_points_dir,      &
       & stress_expansion_order    = 2,                        &
       & stress_subspace_coupling  = stress_subspace_coupling, &
       & vscf_basis_functions_only = .true.,                   &
       & calculation_reader        = calculation_reader,       &
       & logfile                   = logfile                   )
    stress_file = OFile('stress.dat')
    call stress_file%print_lines(stress)
  endif
end subroutine
end module
