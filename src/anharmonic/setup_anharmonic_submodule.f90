submodule (caesar_setup_anharmonic_module) caesar_setup_anharmonic_submodule
  use caesar_anharmonic_module
contains

module procedure setup_anharmonic_mode
  output%mode_name = 'setup_anharmonic'
  output%description = 'Sets up anharmonic calculations. Should be run after &
     &calculate_normal_modes. Several arguments are dependent on the details &
     &of the vibrational potential; this potential can be investigated using &
     &map_potential.'
  output%keywords = [                                                         &
  & KeywordData( 'harmonic_path',                                             &
  &              'harmonic_path is the path to the directory where harmonic &
  &calculations were run.',                                                   &
  &              default_value='.',                                           &
  &              is_path=.true.),                                             &
  & KeywordData( 'q-point_grid',                                              &
  &              'q-point_grid is the number of q-points in each direction in &
  &the Monkhorst-Pack grid on which the anharmonic potential will be &
  &calculated. This should be specified as three integers separated by &
  &spaces.'),                                                                 &
  & KeywordData( 'potential_representation',                                  &
  &              'potential_representation specifies the representation of &
  &the potential which will be used for calculations. Options are &
  &"polynomial".',                                                            &
  &              default_value='polynomial'),                                 &
  & KeywordData( 'max_subspace_coupling',                                     &
  &              'Each term in the potential will contain modes from up to &
  &max_subspace_coupling separate degenerate subspaces. max_subspace_coupling &
  &must be at least 1. This parameter should be converged, although 1 will &
  &often be sufficient.',                                                     &
  &              default_value='1'),                                          &
  & KeywordData( 'max_q-point_coupling',                                      &
  &              'Each term in the potential will contain modes at up to &
  &max_q-point_coupling separate q-points. max_q-point_coupling must be at &
  &least 1. This parameter should be converged, although 1 will often be &
  &sufficient.',                                                              &
  &              default_value='1'),                                          &
  & KeywordData( 'potential_expansion_order',                                 &
  &              'potential_expansion_order is the order up to which the &
  &potential is expanded. e.g. if potential_expansion_order=4 then terms up &
  &to and including u^4 are included. Must be at least 2, and at least as &
  &large as max_subspace_coupling and max_q-point coupling. This parameter &
  &should be converged, although 4 will often be sufficient.',                &
  &              default_value='4'),                                          &
  & KeywordData( 'vscf_basis_functions_only',                                 &
  &              'vscf_basis_functions_only specifies that the potential will &
  &only be expanded in terms of basis functions which are relevant to vscf. &
  &Non-vscf basis functions only affect anharmonic interpolation, and should &
  &usually be safe to neglect without overly affecting accuracy.',            &
  &              default_value='true'),                                       &
  & KeywordData( 'max_energy_of_displacement',                                &
  &              'max_energy_of_displacement, E_max, and &
  &maximum_displacement, u_max, define how far the system is displaced along &
  &each mode when sampling the anharmonic potential. The displacement u along &
  &each mode is bounded by u<u_max and V(u)<E_max, where V is the harmonic &
  &potential, V(u)=0.5*m*(wu)^2, where w is the harmonic frequency of the &
  &mode and m is the effective mass of the mode. E_max and u_max should be &
  &chosen such that the potential is a bounded well, where the boundary &
  &height (the difference between the potential energy at the boundary and &
  &the minimum of the potential) is significantly larger than the maximum &
  &temperature plus the zero-point energy per mode. E_max should be given in &
  &Hartree.',                                                                 &
  &              exclusive_with=[str('frequency_of_max_displacement')]),      &
  & KeywordData( 'maximum_displacement',                                      &
  &              'maximum_displacement, u_max, and &
  &max_energy_of_displacement, E_max, define how far the system is displaced &
  &along each mode when sampling the anharmonic potential. The displacement u &
  &along each mode is bounded by u<u_max and V(u)<E_max, where V is the &
  &harmonic potential, V(u)=0.5*m*(wu)^2, where w is the harmonic frequency &
  &of the mode and m is the effective mass of the mode. E_max and u_max &
  &should be chosen such that the potential is a bounded well, where the &
  &boundary height (the difference between the potential energy at the &
  &boundary and the minimum of the potential) is significantly larger than &
  &the maximum temperature plus the zero-point energy per mode. u_max should &
  &be given in Bohr.'),                                                       &
  & KeywordData( 'frequency_of_max_displacement',                             &
  &              'frequency_of_max_displacement is deprecated.',              &
  &              exclusive_with=[str('max_energy_of_displacement')]),         &
  & KeywordData( 'use_forces',                                                &
  &              'use_forces specifies whether or not each single-point &
  &calculation will produce forces which can be used for fitting the &
  &potential. It is only worth using forces if they are calculated using a &
  &faster method than finite differences.',                                   &
  &              default_value='true'),                                       &
  & KeywordData( 'energy_to_force_ratio',                                     &
  &              'energy_to_force_ratio is the ratio of how penalised &
  &deviations in energy are compared to deviations in forces when the &
  &potential is being fitted. This should be given in units of (Hartree &
  &per primitive cell) divided by (Hartree per bohr). Due to the use of &
  &mass-weighted co-ordinates, in systems containing different elements &
  &forces along modes with higher contributions from heavier elements will &
  &be weighted less than this. Increasing energy_to_force_ratio causes the &
  &fit to focus more on the forces and less on the energies. If the potential &
  &fit is good, changing the value of energy_to_force_ratio should have &
  &little impact on the potential fit.',                                      &
  &              default_value='1'),                                          &
  & KeywordData( 'use_hessians',                                              &
  &              'use_hessians specifies whether or not each single-point &
  &calculation will produce a Hessian (the second derivatives of the &
  &potential) or an equivalent, such as dynamical matrices or phonons, which &
  &can be used for fitting the potential. Only calculations at the Gamma &
  &point of each single point calculation are used. It is only worth using &
  &Hessians if they are calculated using a faster method than finite &
  &differences.',                                                             &
  &              default_value='false'),                                      &
  & KeywordData( 'calculate_stress',                                          &
  &              'calculate_stress specifies whether or not to calculate &
  &stress and pressure. If calculate_stress is true then all electronic &
  &structure calculations must produce a stress tensor or virial tensor. If &
  &calculate_stress is true then potential_expansion_order must be at least &
  &4, otherwise insufficient sampling points are generated to fit the full &
  &stress.',                                                                  &
  &              default_value='true')                                        ]
  output%main_subroutine => setup_anharmonic_subroutine
end procedure

module procedure setup_anharmonic_subroutine
  ! User inputs.
  type(String)          :: harmonic_path
  integer               :: qpoint_grid(3)
  type(String)          :: potential_representation
  integer               :: potential_expansion_order
  integer               :: max_subspace_coupling
  integer               :: max_qpoint_coupling
  logical               :: vscf_basis_functions_only
  real(dp)              :: maximum_displacement
  real(dp), allocatable :: max_energy_of_displacement
  real(dp), allocatable :: frequency_of_max_displacement
  logical               :: use_forces
  real(dp)              :: energy_to_force_ratio
  logical               :: use_hessians
  logical               :: calculate_stress
  
  real(dp) :: weighted_energy_force_ratio
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Harmonic data.
  type(StructureData)                :: structure
  type(StructureData)                :: harmonic_supercell
  type(QpointData),      allocatable :: harmonic_qpoints(:)
  type(DynamicalMatrix), allocatable :: harmonic_dynamical_matrices(:)
  type(ComplexMode),     allocatable :: harmonic_complex_modes(:,:)
  
  ! Data common to all potential representations.
  type(InterpolatedSupercell) :: interpolated_supercell
  type(MaxDisplacement)       :: max_displacement
  type(AnharmonicData)        :: anharmonic_data
  
  ! The potential.
  type(PotentialPointer) :: potential
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
  ! Directories.
  type(String) :: qpoint_dir
  type(String) :: sampling_points_dir
  
  ! Input files.
  type(IFile) :: structure_file
  type(IFile) :: harmonic_supercell_file
  type(IFile) :: harmonic_qpoints_file
  type(IFile) :: harmonic_dynamical_matrices_file
  type(IFile) :: harmonic_complex_modes_file
  
  ! Output files.
  type(OFile) :: logfile
  type(OFile) :: anharmonic_data_file
  type(OFile) :: calculation_directories_file
  
  ! Temporary variables.
  integer :: i,ialloc
  
  ! ----------------------------------------------------------------------
  ! Read in inputs.
  ! ----------------------------------------------------------------------
  
  ! Parse new user inputs.
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split_line(arguments%value('q-point_grid')))
  potential_representation = arguments%value('potential_representation')
  potential_expansion_order = int(arguments%value('potential_expansion_order'))
  max_subspace_coupling = int(arguments%value('max_subspace_coupling'))
  max_qpoint_coupling = int(arguments%value('max_q-point_coupling'))
  vscf_basis_functions_only = &
     & lgcl(arguments%value('vscf_basis_functions_only'))
  maximum_displacement = dble(arguments%value('maximum_displacement'))
  if (arguments%is_set('frequency_of_max_displacement')) then
    frequency_of_max_displacement = &
       & dble(arguments%value('frequency_of_max_displacement'))
  else
    max_energy_of_displacement = &
       & dble(arguments%value('max_energy_of_displacement'))
  endif
  use_forces = lgcl(arguments%value('use_forces'))
  energy_to_force_ratio = dble(arguments%value('energy_to_force_ratio'))
  use_hessians = lgcl(arguments%value('use_hessians'))
  calculate_stress = lgcl(arguments%value('calculate_stress'))
  
  if (potential_expansion_order<2) then
    call print_line(ERROR//': potential_expansion_order must be at least 2.')
    call quit()
  elseif (max_subspace_coupling<1) then
    call print_line(ERROR//': max_subspace_coupling must be at least 1.')
    call quit()
  elseif (max_qpoint_coupling<1) then
    call print_line(ERROR//': max_q-point_coupling must be at least 1.')
    call quit()
  elseif (calculate_stress .and. potential_expansion_order<4) then
    call print_line(ERROR//': calculate_stress may only be true if &
       &potential_expansion_order is at least 4.')
    call quit()
  endif
  
  ! Read setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic_mode())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  call setup_harmonic_arguments%write_file('setup_harmonic.used_settings')
  
  ! ----------------------------------------------------------------------
  ! Read in harmonic data.
  ! ----------------------------------------------------------------------
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  harmonic_supercell_file = IFile(harmonic_path//'/large_supercell.dat')
  harmonic_supercell = StructureData(harmonic_supercell_file%lines())
  
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
  
  ! ----------------------------------------------------------------------
  ! Generate and write out data common to all potential representations.
  ! ----------------------------------------------------------------------
  
  ! Interpolate the supercell, and construct data at the interpolated q-points.
  call print_line('Generating anharmonic supercell.')
  logfile = OFile('setup_anharmonic_logfile.dat')
  interpolated_supercell = InterpolatedSupercell( &
                   & qpoint_grid,                 &
                   & structure,                   &
                   & harmonic_supercell,          &
                   & harmonic_qpoints,            &
                   & harmonic_dynamical_matrices, &
                   & harmonic_complex_modes,      &
                   & logfile                      )
  
  ! Construct the maximum displacement in normal-mode co-ordinates.
  max_displacement = MaxDisplacement(                                 &
     & maximum_displacement          = maximum_displacement,          &
     & structure                     = structure,                     &
     & frequency_of_max_displacement = frequency_of_max_displacement, &
     & max_energy_of_displacement    = max_energy_of_displacement     )
  
  ! Construct and collate common data.
  anharmonic_data = AnharmonicData(                           &
     & structure                 = structure,                 &
     & interpolated_supercell    = interpolated_supercell,    &
     & potential_expansion_order = potential_expansion_order, &
     & max_subspace_coupling     = max_subspace_coupling,     &
     & max_qpoint_coupling       = max_qpoint_coupling,       &
     & vscf_basis_functions_only = vscf_basis_functions_only, &
     & max_displacement          = max_displacement           )
  
  ! Write out anharmonic data.
  call print_line('Writing anharmonic data.')
  anharmonic_data_file = OFile('anharmonic_data.dat')
  call anharmonic_data_file%print_lines(anharmonic_data)
  
  ! ----------------------------------------------------------------------
  ! Generate and write out sampling points.
  ! ----------------------------------------------------------------------
  ! Make a directory for sampling points.
  sampling_points_dir = 'sampling_points'
  call mkdir(sampling_points_dir)
  
  ! Initialise potential to the chosen representation.
  call print_line('Initialising potential.')
  if (potential_representation=='polynomial') then
    potential = PotentialPointer(PolynomialPotential(anharmonic_data))
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
  ! Initialise calculation writer.
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  ! Calculate weighted energy to force ratio.
  weighted_energy_force_ratio = &
     &   energy_to_force_ratio  &
     & * sqrt(maxval(anharmonic_data%structure%atoms%mass()))
  
  ! Generate the sampling points which will be used to map out the anharmonic
  !    Born-Oppenheimer surface in the chosen representation.
  call print_line('Generating sampling points.')
  call potential%generate_sampling_points( anharmonic_data,             &
                                         & use_forces,                  &
                                         & weighted_energy_force_ratio, &
                                         & use_hessians,                &
                                         & calculate_stress,            &
                                         & sampling_points_dir,         &
                                         & calculation_writer,          &
                                         & logfile                      )
  
  ! Write out calculation directories to file.
  calculation_directories_file = OFile('calculation_directories.dat')
  call calculation_directories_file%print_lines( &
     & calculation_writer%directories_written()  )
end procedure
end submodule
