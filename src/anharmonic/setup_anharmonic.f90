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
  &grid used for harmonic calculations.')]
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
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  type(String) :: wd
  type(String) :: harmonic_path
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(String)     :: seedname
  type(String)     :: file_type
  integer          :: qpoint_grid(3)
  
  ! Previously calculated data.
  type(StructureData)            :: structure
  type(QpointData),  allocatable :: harmonic_qpoints(:)
  type(ComplexMode), allocatable :: normal_modes(:)
  
  ! Anharmonic q-points and the corresponding supercell.
  type(IntMatrix)                :: anharmonic_supercell_matrix
  type(StructureData)            :: anharmonic_supercell
  type(QpointData),  allocatable :: qpoints(:)
  
  ! Basis functions.
  integer                        :: no_degenerate_spaces
  integer,           allocatable :: mode_qpoints(:)
  type(ComplexMode), allocatable :: degenerate_modes(:)
  integer,           allocatable :: degenerate_mode_qpoints(:)
  type(Polynomial),  allocatable :: degenerate_mode_bases(:)
  
  ! Temporary variables.
  integer :: i,j,k,l,ialloc
  
  ! Parse inputs.
  wd = arguments%value('working_directory')
  harmonic_path = arguments%value('harmonic_path')
  qpoint_grid = int(split(arguments%value('q-point_grid')))
  
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
  
  ! Generate the basis functions for each degenerate set of modes.
  no_degenerate_spaces = maxval(normal_modes%degeneracy_id)
  do i=1,no_degenerate_spaces
    ! Identify each set of degenerate modes, id by id.
    degenerate_modes = normal_modes(filter(normal_modes%degeneracy_id==i))
    degenerate_mode_qpoints = &
       & mode_qpoints(filter(normal_modes%degeneracy_id==i))
    
    ! Ignore degeneracy ids for which there are no modes.
    ! (missing modes are at harmonic q-points which are not anharmonic.)
    if (size(degenerate_modes)==0) then
      cycle
    endif
    
    call print_line(degenerate_modes%degeneracy_id)
  enddo
  
  ! Combine degenerate sets to get full basis functions.
  ! TODO.
  
  ! Use basis functions to generate sampling points.
  ! TODO
  
  ! Write out sampling points.
  ! TODO
end subroutine
end module
