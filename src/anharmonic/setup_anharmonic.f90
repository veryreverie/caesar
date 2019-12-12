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
  type(String) :: harmonic_path
  integer      :: qpoint_grid(3)
  type(String) :: potential_representation
  integer      :: potential_expansion_order
  integer      :: maximum_coupling_order
  logical      :: vscf_basis_functions_only
  real(dp)     :: maximum_displacement
  real(dp)     :: max_energy_of_displacement
  real(dp)     :: frequency_of_max_displacement
  logical      :: use_forces
  real(dp)     :: energy_to_force_ratio
  logical      :: use_hessians
  logical      :: calculate_stress
  
  real(dp) :: weighted_energy_force_ratio
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  
  ! Previously calculated data.
  type(StructureData)            :: structure
  type(QpointData),  allocatable :: harmonic_qpoints(:)
  type(ComplexMode), allocatable :: qpoint_modes(:)
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
  ! Maximum displacement in mass-weighted co-ordinates.
  real(dp) :: minimum_mass
  real(dp) :: maximum_weighted_displacement
  
  ! Anharmonic q-points and the corresponding supercell.
  type(IntMatrix)                :: anharmonic_supercell_matrix
  type(StructureData)            :: anharmonic_supercell
  type(QpointData),  allocatable :: qpoints(:)
  
  ! Degeneracy data.
  type(DegenerateSubspace), allocatable :: degenerate_subspaces(:)
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  
  ! Coupling data.
  type(SubspaceCoupling), allocatable :: subspace_coupling(:)
  
  ! Anharmonic data container.
  type(AnharmonicData) :: anharmonic_data
  
  ! The potential.
  type(PotentialPointer) :: potential
  
  ! Directories.
  type(String) :: qpoint_dir
  type(String) :: sampling_points_dir
  
  ! Input files.
  type(IFile) :: structure_file
  type(IFile) :: harmonic_qpoints_file
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
  integer :: i,j,ialloc
  
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
  
  ! Read in structure and harmonic q-points.
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  harmonic_qpoints_file = IFile(harmonic_path//'/qpoints.dat')
  harmonic_qpoints = QpointData(harmonic_qpoints_file%sections())
  
  ! Calculate weighted energy to force ratio.
  weighted_energy_force_ratio = &
     &   energy_to_force_ratio  &
     & * sqrt(maxval(structure%atoms%mass()))
  
  ! --------------------------------------------------
  ! Initialise calculation writer.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  ! ----------------------------------------------------------------------
  ! Generate setup data common to all potential representations.
  ! ----------------------------------------------------------------------
  
  ! Calculate the maximum mass-weighted displacement from the maximum
  !    displacement. This corresponds to a mode made entirely from the
  !    lightest element moving up to maximum_displacement.
  ! Also calculate whichever of max_energy_of_displacement and
  !    frequency_of_max_displacement was not given by the user, by
  !    E_max = 0.5 * m_min * (w_max*u_max)^2.
  minimum_mass = minval(structure%atoms%mass())
  maximum_weighted_displacement = maximum_displacement * sqrt(minimum_mass)
  if (arguments%is_set('frequency_of_max_displacement')) then
    max_energy_of_displacement = 0.5_dp * ( frequency_of_max_displacement &
                                        & * maximum_weighted_displacement )**2
  else
    frequency_of_max_displacement = sqrt(2*max_energy_of_displacement) &
                                & / maximum_weighted_displacement
  endif
  
  ! Generate anharmonic q-point grid, and the supercell which has all
  !    anharmonic q-points as G-vectors.
  call print_line('Generating anharmonic supercell.')
  anharmonic_supercell_matrix =                                &
     & mat([ qpoint_grid(1), 0             , 0            ,    &
     &       0             , qpoint_grid(2), 0            ,    &
     &       0             , 0             , qpoint_grid(3) ], &
     & 3,3)
  anharmonic_supercell = construct_supercell( structure,                  &
                                            & anharmonic_supercell_matrix )
  qpoints = generate_qpoints(anharmonic_supercell)
  
  ! Read in harmonic normal modes which correspond to anharmonic q-points,
  !    and record which new q-point each corresponds to.
  allocate(complex_modes(0), stat=ialloc); call err(ialloc)
  do i=1,size(qpoints)
    
    ! Identify the matching harmonic q-point.
    j = first(harmonic_qpoints==qpoints(i),default=0)
    if (j==0) then
      call print_line(ERROR//': anharmonic q-point '//qpoints(i)%qpoint// &
         &' is not also a harmonic q-point.')
      call quit()
    endif
    
    ! Re-label the anharmonic q-point to have the same ID as the matching
    !    harmonic q-point.
    qpoints(i)%id = harmonic_qpoints(j)%id
    qpoints(i)%paired_qpoint_id = harmonic_qpoints(j)%paired_qpoint_id
    
    ! Read in the modes at the harmonic q-point.
    qpoint_dir = &
       & harmonic_path//'/qpoint_'//left_pad(j,str(size(harmonic_qpoints)))
    harmonic_complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    qpoint_modes = ComplexMode(harmonic_complex_modes_file%sections())
    
    complex_modes = [complex_modes, qpoint_modes]
  enddo
  
  complex_modes = complex_modes(filter(.not.complex_modes%translational_mode))
  
  ! Calculate real modes from complex modes.
  real_modes = complex_to_real(complex_modes)
  
  ! Retrieve data on how normal modes are grouped into subspaces
  !    of degenerate modes.
  degenerate_subspaces = process_degeneracies(complex_modes)
  
  ! Generate the symmetry operators in each degenerate subspace.
  allocate( degenerate_symmetries(size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%symmetries)
    degenerate_symmetries(i) = DegenerateSymmetry( structure%symmetries(i), &
                                                 & degenerate_subspaces,    &
                                                 & complex_modes,           &
                                                 & qpoints                  )
  enddo
  
  ! Generate all sets of coupled subspaces, up to maximum_coupling_order.
  call print_line('Generating couplings between subspaces.')
  subspace_coupling = generate_coupled_subspaces( degenerate_subspaces, &
                                                & maximum_coupling_order)
  
  ! Load anharmonic data into container.
  anharmonic_data = AnharmonicData( structure,                     &
                                  & anharmonic_supercell,          &
                                  & qpoints,                       &
                                  & complex_modes,                 &
                                  & real_modes,                    &
                                  & degenerate_subspaces,          &
                                  & degenerate_symmetries,         &
                                  & subspace_coupling,             &
                                  & maximum_coupling_order,        &
                                  & potential_expansion_order,     &
                                  & vscf_basis_functions_only,     &
                                  & maximum_weighted_displacement, &
                                  & frequency_of_max_displacement  )
  
  ! ----------------------------------------------------------------------
  ! Write out setup data common to all potential representations.
  ! ----------------------------------------------------------------------
  call print_line('Writing common data.')
  ! Write out anharmonic supercell and q-points.
  anharmonic_supercell_file = OFile('anharmonic_supercell.dat')
  call anharmonic_supercell_file%print_lines(anharmonic_supercell)
  
  anharmonic_qpoints_file = OFile('anharmonic_qpoints.dat')
  call anharmonic_qpoints_file%print_lines(qpoints, separating_line='')
  
  ! Write out complex and real normal modes.
  complex_modes_file = OFile('complex_modes.dat')
  call complex_modes_file%print_lines(complex_modes, separating_line='')
  
  real_modes_file = OFile('real_modes.dat')
  call real_modes_file%print_lines(real_modes, separating_line='')
  
  ! Write out subspaces and subspace coupling.
  subspaces_file = OFile('degenerate_subspaces.dat')
  call subspaces_file%print_lines(degenerate_subspaces, separating_line='')
  
  coupling_file = OFile('subspace_coupling.dat')
  call coupling_file%print_lines(subspace_coupling)
  
  ! Write out symmetries.
  symmetry_file = OFile('symmetries.dat')
  call symmetry_file%print_lines(degenerate_symmetries, separating_line='')
  
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
    potential = PotentialPointer(                       &
       & PolynomialPotential(potential_expansion_order) )
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
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
