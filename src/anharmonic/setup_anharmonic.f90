! ======================================================================
! Sets up anharmonic calculations.
! ======================================================================
module setup_anharmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
  implicit none
  
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_anharmonic_keywords() result(keywords)
  use keyword_module
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
  & KeywordData( 'potential_expansion_order',                                 &
  &              'potential_expansion_order is the order up to which the &
  &potential is expanded. e.g. if potential_expansion_order=4 then terms up &
  &to and including u^4 are included.'),                                      &
  & KeywordData( 'coupling_order',                                            &
  &              'coupling_order is the maximum number of degenerate &
  &subspaces which are coupled together.') ]
end function

function setup_anharmonic_mode() result(output)
  use caesar_modes_module
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
  use dictionary_module
  use linear_algebra_module
  use setup_harmonic_module
  use structure_module
  use qpoints_module
  use normal_mode_module
  use polynomial_module
  use generate_qpoints_module
  use logic_module
  use integer_arrays_module
  use normal_mode_symmetry_module
  use symmetry_module
  use ofile_module
  use coupling_module
  use basis_function_module
  use degeneracy_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String) :: wd
  type(String) :: harmonic_path
  integer      :: potential_expansion_order
  integer      :: coupling_order
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  integer          :: qpoint_grid(3)
  
  ! Previously calculated data.
  type(StructureData)            :: structure
  type(QpointData),  allocatable :: harmonic_qpoints(:)
  type(ComplexMode), allocatable :: normal_modes(:)
  integer,           allocatable :: mode_qpoints(:)
  
  ! Anharmonic q-points and the corresponding supercell.
  type(IntMatrix)                :: anharmonic_supercell_matrix
  type(StructureData)            :: anharmonic_supercell
  type(QpointData),  allocatable :: qpoints(:)
  
  ! Degeneracy data.
  type(DegenerateModes), allocatable :: degenerate_subspaces(:)
  
  ! Symmetries.
  type(SymmetryOperator)           :: symmetry
  type(ComplexMatrix), allocatable :: symmetries(:,:)
  
  ! Coupling data.
  type(CoupledSubspaces), allocatable :: couplings(:)
  
  ! Basis functions.
  type(Polynomial), allocatable :: basis_functions(:)
  
  ! Logfile.
  type(OFile) :: logfile
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! Parse inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split(arguments%value('q-point_grid')))
  potential_expansion_order = int(arguments%value('potential_expansion_order'))
  coupling_order = int(arguments%value('coupling_order'))
  
  ! Retrieve previous data.
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file( &
     & harmonic_path//'/setup_harmonic.used_settings')
  seedname = setup_harmonic_arguments%value('seedname')
  file_type = setup_harmonic_arguments%value('file_type')
  
  structure = read_structure_file(harmonic_path//'/structure.dat')
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
                                            & calculate_symmetries=.false.)
  call write_structure_file( anharmonic_supercell, &
                           & wd//'/anharmonic_supercell.dat')
  qpoints = generate_qpoints(anharmonic_supercell)
  call write_qpoints_file(qpoints, wd//'/qpoints.dat')
  
  ! Read in normal modes at anharmonic q-points, and record which new q-point
  !    each corresponds to.
  allocate( normal_modes(size(qpoints)*structure%no_modes), &
          & mode_qpoints(size(qpoints)*structure%no_modes), &
          & stat=ialloc); call err(ialloc)
  l = 0
  do i=1,size(qpoints)
    j = first(harmonic_qpoints==qpoints(i))
    
    if (j==0) then
      call print_line(ERROR//': anharmonic q-point '//qpoints(i)%qpoint// &
         &' is not also a harmonic q-point.')
      stop
    endif
    
    do k=1,structure%no_modes
      l = l+1
      normal_modes(l) = ComplexMode(                              &
         & harmonic_path                                       // &
         & '/qpoint_'//left_pad(j,str(size(harmonic_qpoints))) // &
         & '/mode_'//left_pad(k,str(structure%no_modes))//'.dat')
      mode_qpoints(l) = i
    enddo
  enddo
  
  degenerate_subspaces = process_degeneracies(normal_modes,mode_qpoints)
  
  ! Generate the symmetry operators in each subspace.
  allocate( symmetries( size(degenerate_subspaces),  &
          &             size(structure%symmetries)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(structure%symmetries)
    symmetry = structure%symmetries(i)
    do j=1,size(degenerate_subspaces)
      symmetries(j,i) = calculate_symmetry_in_normal_coordinates( &
                   & degenerate_subspaces(j)%modes(normal_modes), &
                   & degenerate_subspaces(j)%qpoints(qpoints),    &
                   & symmetry,                                    &
                   & logfile)
    enddo
  enddo
  
  ! Generate all possible couplings between subspaces.
  couplings = generate_subspace_coupling( degenerate_subspaces,      &
                                        & potential_expansion_order, &
                                        & coupling_order)
  call print_line(size(couplings))
  
  ! Generate basis functions at each coupling.
  basis_functions = [Polynomial::]
  do i=1,size(couplings)
    basis_functions = [ basis_functions,                                &
                    &   generate_basis_functions( couplings(i),         &
                    &                             normal_modes,         &
                    &                             qpoints,              &
                    &                             degenerate_subspaces, &
                    &                             symmetries )          &
                    & ]
  enddo
  
  ! Combine degenerate sets to get full basis functions.
  ! TODO.
  
  ! Use basis functions to generate sampling points.
  ! TODO
  
  ! Write out sampling points.
  ! TODO
end subroutine
end module
