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
  
  public :: setup_anharmonic_keywords
  public :: setup_anharmonic_mode
  public :: setup_anharmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_anharmonic_keywords() result(keywords)
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
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
end function

function setup_anharmonic_mode() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'setup_anharmonic'
  output%description = 'Sets up anharmonic calculations. Should be run after &
     &calculate_harmonic.'
  output%keywords = setup_anharmonic_keywords()
  output%main_subroutine => setup_anharmonic
end function

! ----------------------------------------------------------------------
! The main program.
! ----------------------------------------------------------------------
subroutine setup_anharmonic(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String) :: wd
  type(String) :: harmonic_path
  integer      :: potential_expansion_order
  integer      :: maximum_coupling_order
  logical      :: vscf_basis_functions_only
  real(dp)     :: maximum_displacement
  real(dp)     :: frequency_of_max_displacement
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  integer          :: qpoint_grid(3)
  real(dp)         :: symmetry_precision
  
  ! Previously calculated data.
  type(StructureData)            :: structure
  type(QpointData),  allocatable :: harmonic_qpoints(:)
  type(ComplexMode), allocatable :: complex_modes(:)
  integer,           allocatable :: mode_qpoints(:)
  type(RealMode),    allocatable :: real_modes(:)
  
  ! Anharmonic q-points and the corresponding supercell.
  type(IntMatrix)                :: anharmonic_supercell_matrix
  type(StructureData)            :: anharmonic_supercell
  type(QpointData),  allocatable :: qpoints(:)
  
  ! Degeneracy data.
  type(DegenerateModes), allocatable :: degenerate_subspaces(:)
  
  ! Symmetries.
  type(DegenerateSymmetry), allocatable :: degenerate_symmetries(:)
  
  ! Coupling data.
  type(CoupledSubspaces), allocatable :: coupled_subspaces(:)
  type(SubspaceMonomial), allocatable :: subspace_monomials(:)
  
  ! Basis functions.
  type(BasisFunctions), allocatable :: basis_functions(:)
  
  ! Sampling points.
  real(dp)                          :: maximum_weighted_displacement
  type(SamplingPoints), allocatable :: sampling_points(:)
  
  ! Supercell data.
  type(IntMatrix)     :: supercell_matrix
  type(StructureData) :: supercell
  
  ! Displacement data.
  type(VscfRvectors), allocatable :: vscf_rvectors(:)
  type(CartesianDisplacement)     :: displacement
  type(StructureData)             :: displaced_structure
  
  ! Potential data.
  type(PotentialPointer) :: potential
  
  ! Directories and files.
  type(String) :: qpoint_dir
  type(String) :: coupling_dir
  type(String) :: sampling_points_dir
  type(String) :: sampling_dir
  type(String) :: vscf_rvector_dir
  type(String) :: input_filename
  
  ! Input files.
  type(IFile)                    :: harmonic_qpoints_file
  type(IFile)                    :: harmonic_complex_modes_file
  type(StringArray), allocatable :: file_sections(:)
  
  ! Output files.
  type(OFile) :: logfile
  type(OFile) :: anharmonic_qpoints_file
  type(OFile) :: complex_modes_file
  type(OFile) :: real_modes_file
  type(OFile) :: coupling_file
  type(OFile) :: basis_function_file
  type(OFile) :: sampling_points_file
  type(OFile) :: vscf_rvectors_file
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! ----------------------------------------------------------------------
  ! Read in inputs.
  ! ----------------------------------------------------------------------
  
  ! Parse new user inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split_line(arguments%value('q-point_grid')))
  potential_expansion_order = int(arguments%value('potential_expansion_order'))
  maximum_coupling_order = int(arguments%value('maximum_coupling_order'))
  vscf_basis_functions_only = &
     & lgcl(arguments%value('vscf_basis_functions_only'))
  maximum_displacement = dble(arguments%value('maximum_displacement'))
  frequency_of_max_displacement = &
     & dble(arguments%value('frequency_of_max_displacement'))
  
  ! Retrieve data from previous stages.
  setup_harmonic_arguments = Dictionary(setup_harmonic_keywords())
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & symmetry_precision)
  
  harmonic_qpoints_file = IFile(harmonic_path//'/qpoints.dat')
  file_sections = split_into_sections(harmonic_qpoints_file%lines())
  allocate( harmonic_qpoints(size(file_sections)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(harmonic_qpoints)
    harmonic_qpoints(i) = file_sections(i)
  enddo
  
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
  anharmonic_supercell = construct_supercell( structure,                   &
                                            & anharmonic_supercell_matrix, &
                                            & symmetry_precision,          &
                                            & calculate_symmetry=.false.)
  qpoints = generate_qpoints(anharmonic_supercell)
  
  ! Read in harmonic normal modes which correspond to anharmonic q-points,
  !    and record which new q-point each corresponds to.
  allocate( complex_modes(size(qpoints)*structure%no_modes), &
          & real_modes(size(qpoints)*structure%no_modes),    &
          & mode_qpoints(size(qpoints)*structure%no_modes),  &
          & stat=ialloc); call err(ialloc)
  l = 0
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
    file_sections = split_into_sections(harmonic_complex_modes_file%lines())
    do k=1,structure%no_modes
      l = l+1
      complex_modes(l) = file_sections(k)
      mode_qpoints(l) = i
    enddo
  enddo
  
  ! Calculate real modes from complex modes.
  real_modes = complex_to_real(complex_modes)
  
  ! Retrieve data on how normal modes are grouped into subspaces
  !    of degenerate modes.
  degenerate_subspaces = process_degeneracies(complex_modes,mode_qpoints)
  
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
  coupled_subspaces = generate_coupled_subspaces( degenerate_subspaces, &
                                                & maximum_coupling_order)
  
  ! ----------------------------------------------------------------------
  ! Write out setup data common to all potential representations.
  ! ----------------------------------------------------------------------
  
  ! Write out complex and real normal modes.
  complex_modes_file = OFile(wd//'/complex_modes.dat')
  call complex_modes_file%print_lines(complex_modes,separating_line='')
  
  real_modes_file = OFile(wd//'/real_modes.dat')
  call real_modes_file%print_lines(real_modes,separating_line='')
  
  ! Write out anharmonic supercell and q-points.
  call write_structure_file( anharmonic_supercell, &
                           & wd//'/anharmonic_supercell.dat')
  
  anharmonic_qpoints_file = OFile(wd//'/qpoints.dat')
  call anharmonic_qpoints_file%print_lines(qpoints,separating_line='')
  
  ! Write out coupling and basis functions.
  coupling_file = OFile(wd//'/coupling.dat')
  call coupling_file%print_lines(coupled_subspaces)
  
  ! ----------------------------------------------------------------------
  ! Generate and write out sampling points.
  ! ----------------------------------------------------------------------
  ! Make a directory for sampling points.
  sampling_points_dir = wd//'/sampling_points'
  call mkdir(sampling_points_dir)
  
  ! Initialise potential to the chosen representation.
  potential = PolynomialPotential(potential_expansion_order)
  
  ! Generate the sampling points which will be used to map out the anharmonic
  !    Born-Oppenheimer surface in the chosen representation.
  call potential%generate_sampling_points( &
          & sampling_points_dir,           &
          & structure,                     &
          & symmetry_precision,            &
          & complex_modes,                 &
          & real_modes,                    &
          & qpoints,                       &
          & degenerate_subspaces,          &
          & degenerate_symmetries,         &
          & coupled_subspaces,             &
          & vscf_basis_functions_only,     &
          & maximum_weighted_displacement, &
          & frequency_of_max_displacement, &
          & logfile,                       &
          & write_structure_file_lambda)
contains
! Lambda of type WriteLambda to write a structure to file.
! Captures:
!    - file_type
!    - seedname
!    - wd
subroutine write_structure_file_lambda(structure,directory)
  implicit none
  
  type(StructureData), intent(in) :: structure
  type(String),        intent(in) :: directory
  
  type(String) :: input_filename
  
  input_filename = make_input_filename(file_type,seedname)
  call write_structure_file(structure,directory//'/structure.dat')
  call StructureData_to_input_file( &
         & file_type,               &
         & structure,               &
         & wd//'/'//input_filename, &
         & directory//'/'//input_filename)
end subroutine
end subroutine
end module
