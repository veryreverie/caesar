! ======================================================================
! The first stage of Caesar.
! Generates supercells, and prepares harmonic DFT calculations.
! ======================================================================
module setup_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_harmonic_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'file_type',                                                 &
  &              'file_type is the file type which will be used for &
  &single-point energy calculations. Settings are: castep quip.',             &
  &              default_value='castep'),                                     &
  & KeywordData( 'seedname',                                                  &
  &              'seedname is the seedname from which file names are &
  &constructed.'),                                                            &
  & KeywordData( 'q-point_grid',                                              &
  &              'q-point_grid is the number of q-points in each direction &
  &in a Monkhorst-Pack grid. This should be specified as three integers &
  &separated by spaces.'),                                                    &
  & KeywordData( 'symmetry_precision',                                        &
  &              'symmetry_precision is the tolerance at which symmetries &
  &are calculated.',                                                          &
  &              default_value='0.1'),                                        &
  & KeywordData( 'harmonic_displacement',                                     &
  &              'harmonic_displacement is the distance in bohr by which &
  &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
  &surface.',                                                                 &
  &              default_value='0.01') ]
end function

function setup_harmonic_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'setup_harmonic'
  output%description = 'Sets up harmonic calculation. Generates supercells, &
     &and prepares DFT inputs.'
  output%keywords = setup_harmonic_keywords()
  output%main_subroutine => setup_harmonic
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine setup_harmonic(arguments)
  use utils_module, only : mkdir
  use ofile_module
  use linear_algebra_module
  use structure_module
  use group_module
  use qpoints_module
  use dictionary_module
  use input_file_module
  use generate_supercells_module
  use generate_qpoints_module
  use unique_directions_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User input variables.
  type(String) :: wd
  type(String) :: file_type
  type(String) :: seedname
  
  ! User input data.
  type(StructureData) :: structure
  integer             :: grid(3)
  real(dp)            :: symmetry_precision
  real(dp)            :: harmonic_displacement
  
  ! Supercell data.
  type(IntMatrix)                  :: large_supercell_matrix
  type(StructureData)              :: large_supercell
  type(QpointData),    allocatable :: qpoints(:)
  type(StructureData), allocatable :: supercells(:)
  integer                          :: no_supercells
  type(StructureData)              :: supercell
  
  ! Directories.
  type(String) :: sdir
  type(String) :: path
  
  ! Perturbation direction information.
  type(UniqueDirection), allocatable :: unique_directions(:)
  integer                            :: atom
  type(String)                       :: atom_string
  type(String)                       :: direction
  type(RealVector)                   :: displacement
  type(RealVector)                   :: position
  
  ! Temporary variables.
  integer :: i,j
  
  ! Files.
  type(String) :: input_filename
  type(OFile)  :: no_supercells_file
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  grid = int(split(arguments%value('q-point_grid')))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  
  ! Check dft code is supported
  if (file_type/='castep' .and. file_type/='quip') then
    call print_line('')
    call print_line('Error: The file type '//file_type//' is not currently &
       & supported.')
    call print_line('Please choose one of: castep quip.')
    stop
  endif
  
  ! Check dft input files exists.
  input_filename = make_input_filename(file_type,seedname)
  input_filename = wd//'/'//input_filename
  if (.not. file_exists(input_filename)) then
    call print_line('Error: The input file '//input_filename// &
       &' does not exist.')
    stop
  endif
  
  ! --------------------------------------------------
  ! Read in input files.
  ! --------------------------------------------------
  structure = input_file_to_StructureData( file_type,      &
                                         & input_filename, &
                                         & symmetry_precision)
  call write_structure_file(structure,wd//'/structure.dat')
  
  ! --------------------------------------------------
  ! Generate large supercell, for which all q-points are G-vectors.
  ! --------------------------------------------------
  large_supercell_matrix = mat([ grid(1), 0      , 0      , &
                               & 0      , grid(2), 0      , &
                               & 0      , 0      , grid(3)  ], 3,3)
  large_supercell = construct_supercell( structure,              &
                                       & large_supercell_matrix, &
                                       & calculate_symmetries=.false.)
  call write_structure_file(large_supercell, wd//'/large_supercell.dat')
  
  ! --------------------------------------------------
  ! Generate supercells.
  ! --------------------------------------------------
  ! Generate q-points in Monkhorst-Pack grid.
  qpoints = generate_qpoints(large_supercell)
  call write_qpoints_file(qpoints, wd//'/qpoints.dat')
  
  ! Generate non-diagonal supercells.
  supercells = generate_supercells(structure, qpoints, symmetry_precision)
  no_supercells = size(supercells)
  no_supercells_file = wd//'/no_supercells.dat'
  call no_supercells_file%print_line(no_supercells)
  
  ! Loop over supercells.
  do i=1,no_supercells
    sdir=wd//'/Supercell_'//left_pad(i,str(no_supercells))
    
    call mkdir(sdir)
    
    ! Write out structure files.
    supercell = supercells(i)
    call write_structure_file(supercell, sdir//'/structure.dat')
    
    ! Calculate which forces need calculating.
    unique_directions = calculate_unique_directions( supercell, &
                                                   & harmonic_displacement)
    call write_unique_directions_file( unique_directions, &
                                     & sdir//'/unique_directions.dat')
    
    ! --------------------------------------------------
    ! Write energy and force calculation input files.
    ! --------------------------------------------------
    do j=1,size(unique_directions)
      atom = unique_directions(j)%atom_id
      direction = unique_directions(j)%direction
      displacement = unique_directions(j)%displacement
      
      atom_string = left_pad(atom, str(structure%no_atoms))
      path = sdir//'/atom.'//atom_string//'.'//direction
      
      ! Make harmonic run directories.
      call mkdir(path)
      
      ! Move relevant atom.
      position = supercell%atoms(atom)%cartesian_position()
      call supercell%atoms(atom)%set_cartesian_position( position &
                                                     & + displacement)
        
      ! Write energy and force calculation input file.
      input_filename = make_input_filename(file_type,seedname)
      call StructureData_to_input_file(     &
                 & file_type,               &
                 & supercell,               &
                 & wd//'/'//input_filename, &
                 & path//'/'//input_filename)
      
      ! Reset moved atom.
      call supercell%atoms(atom)%set_cartesian_position(position)
    enddo
  enddo
end subroutine
end module
