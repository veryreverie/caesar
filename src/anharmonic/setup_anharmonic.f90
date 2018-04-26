! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use common_module
  
  use setup_harmonic_module
  
  use degeneracy_module
  use coupled_subspaces_module
  use subspace_monomial_module
  use mode_monomial_module
  use degenerate_symmetry_module
  use basis_function_module
  use basis_functions_module
  use sampling_points_module
  implicit none
  
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
  &atoms will be displaced less far than this.') ]
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
  
  ! Sampling points and displacement data.
  real(dp)                          :: maximum_weighted_displacement
  type(SamplingPoints), allocatable :: sampling_points(:)
  
  ! Directories and files.
  type(OFile)               :: coupling_file
  type(String)              :: max_subspace_id
  type(String), allocatable :: coupling_strings(:)
  type(String)              :: coupling_dir
  
  ! Logfile.
  type(OFile) :: logfile
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! Parse inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split(arguments%value('q-point_grid')))
  potential_expansion_order = int(arguments%value('potential_expansion_order'))
  maximum_coupling_order = int(arguments%value('maximum_coupling_order'))
  vscf_basis_functions_only = &
     & lgcl(arguments%value('vscf_basis_functions_only'))
  maximum_displacement = dble(arguments%value('maximum_displacement'))
  
  ! Retrieve previous data.
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  symmetry_precision = &
     & dble(setup_harmonic_arguments%value('symmetry_precision'))
  
  structure = read_structure_file( harmonic_path//'/structure.dat', &
                                 & symmetry_precision)
  harmonic_qpoints = read_qpoints_file(harmonic_path//'/qpoints.dat')
  
  ! Open logifile.
  logfile = wd//'/setup_anharmonic_logfile.dat'
  
  ! Generate new q-point grid, and the supercell which has all anharmonic
  !    q-points as G-vectors.
  anharmonic_supercell_matrix =                                &
     & mat([ qpoint_grid(1), 0             , 0            ,    &
     &       0             , qpoint_grid(2), 0            ,    &
     &       0             , 0             , qpoint_grid(3) ], &
     & 3,3)
  anharmonic_supercell = construct_supercell( structure,                   &
                                            & anharmonic_supercell_matrix, &
                                            & symmetry_precision,          &
                                            & calculate_symmetry=.false.)
  call write_structure_file( anharmonic_supercell, &
                           & wd//'/anharmonic_supercell.dat')
  qpoints = generate_qpoints(anharmonic_supercell)
  call write_qpoints_file(qpoints, wd//'/qpoints.dat')
  
  ! Read in normal modes at anharmonic q-points, and record which new q-point
  !    each corresponds to.
  allocate( complex_modes(size(qpoints)*structure%no_modes), &
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
    
    do k=1,structure%no_modes
      l = l+1
      complex_modes(l) = ComplexMode(                             &
         & harmonic_path                                       // &
         & '/qpoint_'//left_pad(j,str(size(harmonic_qpoints))) // &
         & '/mode_'//left_pad(k,str(structure%no_modes))//'.dat')
      mode_qpoints(l) = i
    enddo
  enddo
  
  degenerate_subspaces = process_degeneracies(complex_modes,mode_qpoints)
  max_subspace_id = str(maxval(degenerate_subspaces%id))
  
  ! Generate the symmetry operators in each subspace.
  allocate( degenerate_symmetries(size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%symmetries)
    degenerate_symmetries(i) = DegenerateSymmetry( structure%symmetries(i), &
                                                 & degenerate_subspaces,    &
                                                 & complex_modes,           &
                                                 & qpoints,                 &
                                                 & logfile)
  enddo
  
  ! Construct real modes from complex modes.
  real_modes = complex_to_real(complex_modes)
  
  ! Calculate the maximum mass-weighted displacement from the maximum
  !    displacement. This corresponds to a mode made entirely from the
  !    lightest element moving up to maximum_displacement.
  maximum_weighted_displacement = maximum_displacement &
                              & * sqrt(minval(structure%atoms%mass()))
  
  ! Generate all sets of coupled subspaces, up to maximum_coupling_order.
  coupled_subspaces = generate_coupled_subspaces( degenerate_subspaces, &
                                                & maximum_coupling_order)
  
  ! Loop over subspace couplings, generating basis functions and sampling
  !    points for each.
  allocate( basis_functions(size(coupled_subspaces)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(coupled_subspaces)
    ! Generate the set of subspace monomials corresponding to the subspace
    !    coupling.
    ! e.g. the coupling [1,2] might have monomials [1,2], [1,1,2] and [1,2,2].
    subspace_monomials = generate_subspace_monomials( &
                              & coupled_subspaces(i), &
                              & potential_expansion_order)
    ! Generate basis functions at each coupling.
    basis_functions(i) = generate_basis_functions( subspace_monomials,        &
                                                 & structure,                 &
                                                 & complex_modes,             &
                                                 & real_modes,                &
                                                 & qpoints,                   &
                                                 & degenerate_subspaces,      &
                                                 & degenerate_symmetries,     &
                                                 & vscf_basis_functions_only, &
                                                 & logfile)
    ! Generate a set of sampling points from which the basis functions can
    !    be constructed.
    sampling_points(i) = generate_sampling_points( &
                   & basis_functions(i)%functions, &
                   & potential_expansion_order,    &
                   & maximum_weighted_displacement)
  enddo
  
  ! Write out coupling and basis functions.
  coupling_file = wd//'/coupling.dat'
  do i=1,size(coupled_subspaces)
    ! Add entry to coupling file for this coupling.
    call coupling_file%print_line(coupled_subspaces(i))
    
    ! Construct directory name from coupled subspace ids.
    coupling_strings = left_pad( coupled_subspaces(i)%ids, &
                               & max_subspace_id)
    coupling_dir = wd//'/coupling_'//join(coupling_strings,delimiter='-')
    
    call mkdir(coupling_dir)
    
    ! Write basis functions to file.
    call basis_functions(i)%write_file(coupling_dir//'/basis_functions.dat')
  enddo
  
  ! Write out sampling points.
  ! TODO
end subroutine
end module
