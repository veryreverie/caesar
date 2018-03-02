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
  integer,          allocatable :: degeneracy_ids(:)
  type(IntArray1D), allocatable :: degenerate_mode_ids(:)
  type(IntArray1D), allocatable :: degenerate_qpoint_ids(:)
  
  ! Coupling data.
  type(IntArray1D), allocatable :: couplings(:)
  
  ! Basis functions.
  type(ComplexMode), allocatable :: degenerate_modes(:)
  integer,           allocatable :: degenerate_mode_qpoints(:)
  type(Polynomial),  allocatable :: degenerate_mode_bases(:)
  
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
  
  ! Identify the modes and corresponding q-points in each subspace.
  degeneracy_ids = normal_modes%degeneracy_id
  degeneracy_ids = degeneracy_ids(set(degeneracy_ids))
  allocate( degenerate_mode_ids(size(degeneracy_ids)), &
          & degenerate_qpoint_ids(size(degeneracy_ids)), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(degeneracy_ids)
    degenerate_modes = normal_modes(filter(normal_modes%degeneracy_id==degeneracy_ids(i)))
    degenerate_mode_ids(i) = degenerate_modes%id
    degenerate_qpoint_ids(i) = mode_qpoints(filter(normal_modes%degeneracy_id==degeneracy_ids(i)))
  enddo
  
  ! Generate all possible couplings between subspaces.
  allocate( couplings(size(degeneracy_ids)**potential_expansion_order), &
          & stat=ialloc); call err(ialloc)
  do i=1,size(couplings)
    allocate( couplings(i)%i(potential_expansion_order), &
            & stat=ialloc); call err(ialloc)
    do j=1,potential_expansion_order
      k = (i-1) / (size(degeneracy_ids)**(j-1))
      k = modulo(k,size(degeneracy_ids))+1
      couplings(i)%i(j) = degeneracy_ids(k)
    enddo
    call print_line(couplings(i))
  enddo
  
  ! Remove couplings with more than coupling_order coupled subspaces.
  couplings = couplings(filter(couplings,allowed_coupling))
  
  do i=1,size(couplings)
    call print_line(couplings(i))
  enddo
  
  ! Combine degenerate sets to get full basis functions.
  ! TODO.
  
  ! Use basis functions to generate sampling points.
  ! TODO
  
  ! Write out sampling points.
  ! TODO
contains
  ! Lambda for checking if a coupling has coupling_order subspaces or fewer.
  ! Also checks that the coupling is in the sorted order.
  ! Captures coupling_order.
  function allowed_coupling(input) result(output)
    implicit none
    
    class(*), intent(in) :: input
    logical              :: output
    
    select type(input); type is(IntArray1D)
      output = size(set(input%i)) <= coupling_order .and. &
            & all(input%i(sort(input%i)) == input%i)
    end select
  end function
end subroutine
end module
