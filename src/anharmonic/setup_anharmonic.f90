! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use common_module
  
  use anharmonic_common_module
  use potentials_module
  implicit none
  
  private
  
  public :: startup_setup_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_setup_anharmonic()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'setup_anharmonic'
  mode%description = 'Sets up anharmonic calculations. Should be run after &
     &calculate_normal_modes.'
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
  & KeywordData( 'potential_representation',                                  &
  &              'potential_representation specifies the representation of &
  &the potential which will be used for calculations. Options are &
  &"polynomial".',                                                            &
  &              default_value='polynomial'),                                 &
  & KeywordData( 'maximum_coupling_order',                                    &
  &              'maximum_coupling_order is the maximum number of degenerate &
  &subspaces which may be coupled together. Must be at least 1.'),            &
  & KeywordData( 'potential_expansion_order',                                 &
  &              'potential_expansion_order is the order up to which the &
  &potential is expanded. e.g. if potential_expansion_order=4 then terms up &
  &to and including u^4 are included. Must be at least 2, and at least as &
  &large as maximum_coupling_order.'),                                        &
  & KeywordData( 'vscf_basis_functions_only',                                 &
  &              'vscf_basis_functions_only specifies that the potential will &
  &only be expanded in terms of basis functions which are relevant to vscf.', &
  &              default_value='true'),                                       &
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
  &be weighted less than this.'),                                             &
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
  &structure calculations must produce a stress tensor or virial tensor.',    &
  &              default_value='true')                                        ]
  mode%main_subroutine => setup_anharmonic_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_anharmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String)          :: harmonic_path
  integer               :: qpoint_grid(3)
  type(String)          :: potential_representation
  integer               :: potential_expansion_order
  integer               :: maximum_coupling_order
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
  type(IFile) :: harmonic_qpoints_file
  type(IFile) :: harmonic_dynamical_matrices_file
  type(IFile) :: harmonic_complex_modes_file
  
  ! Output files.
  type(OFile) :: logfile
  type(OFile) :: anharmonic_data_file
  type(Ofile) :: anharmonic_supercell_file
  type(OFile) :: anharmonic_qpoints_file
  type(OFile) :: complex_modes_file
  type(OFile) :: real_modes_file
  type(OFile) :: subspaces_file
  type(OFile) :: coupling_file
  type(OFile) :: symmetry_file
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
  maximum_coupling_order = int(arguments%value('maximum_coupling_order'))
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
  
  ! Read setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(CaesarMode('setup_harmonic'))
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
  interpolated_supercell = InterpolatedSupercell( &
                   & qpoint_grid,                 &
                   & structure,                   &
                   & harmonic_qpoints,            &
                   & harmonic_dynamical_matrices, &
                   & harmonic_complex_modes       )
  
  ! Construct the maximum displacement in normal-mode co-ordinates.
  max_displacement = MaxDisplacement(                                 &
     & maximum_displacement          = maximum_displacement,          &
     & structure                     = structure,                     &
     & frequency_of_max_displacement = frequency_of_max_displacement, &
     & max_energy_of_displacement    = max_energy_of_displacement     )
  
  ! Construct and collate common data.
  anharmonic_data = AnharmonicData(                                   &
     & structure                     = structure,                     &
     & interpolated_supercell        = interpolated_supercell,        &
     & max_displacement              = max_displacement,              &
     & potential_expansion_order     = potential_expansion_order,     &
     & maximum_coupling_order        = maximum_coupling_order,        &
     & vscf_basis_functions_only     = vscf_basis_functions_only,     &
     & energy_to_force_ratio         = energy_to_force_ratio          )
  
  call print_line('Writing common data.')
  ! Write out anharmonic supercell and q-points.
  anharmonic_supercell_file = OFile('anharmonic_supercell.dat')
  call anharmonic_supercell_file%print_lines( &
     & anharmonic_data%anharmonic_supercell )
  
  anharmonic_qpoints_file = OFile('anharmonic_qpoints.dat')
  call anharmonic_qpoints_file%print_lines( anharmonic_data%qpoints, &
                                          & separating_line=''       )
  
  ! Write out complex and real normal modes.
  complex_modes_file = OFile('complex_modes.dat')
  call complex_modes_file%print_lines( anharmonic_data%complex_modes, &
                                     & separating_line=''             )
  
  real_modes_file = OFile('real_modes.dat')
  call real_modes_file%print_lines( anharmonic_data%real_modes, &
                                  & separating_line=''          )
  
  ! Write out subspaces and subspace coupling.
  subspaces_file = OFile('degenerate_subspaces.dat')
  call subspaces_file%print_lines( anharmonic_data%degenerate_subspaces, &
                                 & separating_line=''                    )
  
  coupling_file = OFile('subspace_coupling.dat')
  call coupling_file%print_lines(anharmonic_data%subspace_couplings)
  
  ! Write out symmetries.
  symmetry_file = OFile('symmetries.dat')
  call symmetry_file%print_lines( anharmonic_data%degenerate_symmetries, &
                                & separating_line=''                     )
  
  ! Write out anharmonic data.
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
  logfile = OFile('setup_anharmonic_logfile.dat')
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
end subroutine
end module
