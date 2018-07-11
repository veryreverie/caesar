! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use anharmonic_common_module
  use polynomial_module
  implicit none
  
  private
  
  public :: setup_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_anharmonic() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'setup_anharmonic'
  output%description = 'Sets up anharmonic calculations. Should be run after &
     &calculate_harmonic.'
  output%keywords = [                                                         &
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
  &              'maximum_displacement is the largest distance any sampling &
  &point will be from the equilibrium position. maximum_displacement should &
  &be given in Bohr. Due to the use of mass-reduced co-ordinates, in systems &
  &containing different elements modes with higher contributions from heavier &
  &atoms will be displaced less far than this.'),                             &
  & KeywordData( 'frequency_of_max_displacement',                             &
  &              'frequency_of_max_displacement is the frequency, w_min, at &
  &which maximum displacement happens. Displacement along modes with w>w_min &
  &is scaled by sqrt(w_min/w), and displacement along modes with w<w_min &
  &is unscaled. w_min should be given in Hartree.') ]
  output%main_subroutine => setup_anharmonic_subroutine
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_anharmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String) :: wd
  type(String) :: harmonic_path
  integer      :: qpoint_grid(3)
  type(String) :: potential_representation
  integer      :: potential_expansion_order
  integer      :: maximum_coupling_order
  logical      :: vscf_basis_functions_only
  real(dp)     :: maximum_displacement
  real(dp)     :: frequency_of_max_displacement
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  real(dp)         :: symmetry_precision
  
  ! Previously calculated data.
  type(StructureData)            :: structure
  type(QpointData),  allocatable :: harmonic_qpoints(:)
  type(ComplexMode), allocatable :: qpoint_modes(:)
  type(ComplexMode), allocatable :: complex_modes(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
  ! Maximum displacement in mass-weighted co-ordinates.
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
  type(Ofile) :: anharmonic_supercell_file
  type(OFile) :: anharmonic_qpoints_file
  type(OFile) :: complex_modes_file
  type(OFile) :: real_modes_file
  type(OFile) :: subspaces_file
  type(OFile) :: coupling_file
  type(OFile) :: calculation_directories_file
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! ----------------------------------------------------------------------
  ! Read in inputs.
  ! ----------------------------------------------------------------------
  
  ! Parse new user inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split_line(arguments%value('q-point_grid')))
  potential_representation = arguments%value('potential_representation')
  potential_expansion_order = int(arguments%value('potential_expansion_order'))
  maximum_coupling_order = int(arguments%value('maximum_coupling_order'))
  vscf_basis_functions_only = &
     & lgcl(arguments%value('vscf_basis_functions_only'))
  maximum_displacement = dble(arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(arguments%value('frequency_of_max_displacement'))
  
  ! Read setup_harmonic arguments.
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  ! Read in structure and harmonic q-points.
  structure_file = IFile(harmonic_path//'/structure.dat')
  structure = StructureData(structure_file%lines())
  
  harmonic_qpoints_file = IFile(harmonic_path//'/qpoints.dat')
  harmonic_qpoints = QpointData(harmonic_qpoints_file%sections())
  
  ! --------------------------------------------------
  ! Initialise calculation writer.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( working_directory = wd,        &
                                        & file_type         = file_type, &
                                        & seedname          = seedname   )
  
  ! ----------------------------------------------------------------------
  ! Generate setup data common to all potential representations.
  ! ----------------------------------------------------------------------
  
  ! Open logfile.
  logfile = OFile(wd//'/setup_anharmonic_logfile.dat')
  
  ! Calculate the maximum mass-weighted displacement from the maximum
  !    displacement. This corresponds to a mode made entirely from the
  !    lightest element moving up to maximum_displacement.
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! Generate anharmonic q-point grid, and the supercell which has all
  !    anharmonic q-points as G-vectors.
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
  complex_modes = [ComplexMode::]
  do i=1,size(qpoints)
    j = first(harmonic_qpoints==qpoints(i),default=0)
    
    if (j==0) then
      call print_line(ERROR//': anharmonic q-point '//qpoints(i)%qpoint// &
         &' is not also a harmonic q-point.')
      stop
    endif
    
    qpoint_dir = &
       & harmonic_path//'/qpoint_'//left_pad(j,str(size(harmonic_qpoints)))
    harmonic_complex_modes_file = IFile(qpoint_dir//'/complex_modes.dat')
    qpoint_modes = ComplexMode(harmonic_complex_modes_file%sections())
    qpoint_modes%qpoint_id = qpoints(j)%id
    complex_modes = [complex_modes, qpoint_modes]
  enddo
  
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
                                                 & qpoints,                 &
                                                 & logfile)
  enddo
  
  ! Generate all sets of coupled subspaces, up to maximum_coupling_order.
  subspace_coupling = generate_coupled_subspaces( degenerate_subspaces, &
                                                & maximum_coupling_order)
  
  ! ----------------------------------------------------------------------
  ! Write out setup data common to all potential representations.
  ! ----------------------------------------------------------------------
  
  ! Write out anharmonic supercell and q-points.
  anharmonic_supercell_file = OFile(wd//'/anharmonic_supercell.dat')
  call anharmonic_supercell_file%print_lines(anharmonic_supercell)
  
  anharmonic_qpoints_file = OFile(wd//'/qpoints.dat')
  call anharmonic_qpoints_file%print_lines(qpoints,separating_line='')
  
  ! Write out complex and real normal modes.
  complex_modes_file = OFile(wd//'/complex_modes.dat')
  call complex_modes_file%print_lines(complex_modes,separating_line='')
  
  real_modes_file = OFile(wd//'/real_modes.dat')
  call real_modes_file%print_lines(real_modes,separating_line='')
  
  ! Write out subspaces and subspace coupling.
  subspaces_file = OFile(wd//'/degenerate_subspaces.dat')
  call subspaces_file%print_lines(degenerate_subspaces,separating_line='')
  
  coupling_file = OFile(wd//'/subspace_coupling.dat')
  call coupling_file%print_lines(subspace_coupling)
  
  ! ----------------------------------------------------------------------
  ! Generate and write out sampling points.
  ! ----------------------------------------------------------------------
  ! Make a directory for sampling points.
  sampling_points_dir = wd//'/sampling_points'
  call mkdir(sampling_points_dir)
  
  ! Initialise potential to the chosen representation.
  if (potential_representation=='polynomial') then
    potential = PolynomialPotential(potential_expansion_order)
  else
    call print_line( ERROR//': Unrecognised potential representation: '// &
                   & potential_representation)
    call err()
  endif
  
  ! Load anharmonic data into container.
  anharmonic_data = AnharmonicData( structure,                     &
                                  & symmetry_precision,            &
                                  & anharmonic_supercell,          &
                                  & qpoints,                       &
                                  & complex_modes,                 &
                                  & real_modes,                    &
                                  & degenerate_subspaces,          &
                                  & degenerate_symmetries,         &
                                  & subspace_coupling,             &
                                  & vscf_basis_functions_only,     &
                                  & maximum_weighted_displacement, &
                                  & frequency_of_max_displacement )
  
  ! Generate the sampling points which will be used to map out the anharmonic
  !    Born-Oppenheimer surface in the chosen representation.
  call potential%generate_sampling_points( &
                    & anharmonic_data,     &
                    & sampling_points_dir, &
                    & calculation_writer,  &
                    & logfile              )
  
  ! Write out calculation directories to file.
  calculation_directories_file = OFile(wd//'/calculation_directories.dat')
  call calculation_directories_file%print_lines( &
     & calculation_writer%directories_written())
end subroutine
end module
