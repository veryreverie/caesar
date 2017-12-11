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
  &              default_value='0.1')]
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
  
  ! Supercell data.
  type(IntMatrix)                  :: large_supercell_matrix
  type(StructureData)              :: large_supercell
  type(QpointData),    allocatable :: qpoints(:)
  type(StructureData), allocatable :: supercells(:)
  integer                          :: no_supercells
  type(StructureData)              :: supercell
  
  ! Directories.
  type(String) :: sdir
  type(String) :: paths(2)
  
  ! Perturbation direction information.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction_char
  type(RealVector)       :: displacement
  type(String)           :: atom_string
  
  ! Temporary variables.
  integer        :: i,j,k
  
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
  ! Generate q-points in IBZ and non-diagonal supercells.
  qpoints = generate_qpoints(structure, large_supercell)
  
  supercells = generate_supercells(structure, qpoints, symmetry_precision)
  
  call write_qpoints_file(qpoints, wd//'/qpoints.dat')
  
  ! Write no_supercells to file
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
    unique_directions = calculate_unique_directions(supercell)
    call write_unique_directions_file( unique_directions, &
                                     & sdir//'/unique_directions.dat')
    
    ! --------------------------------------------------
    ! Write DFT code input files.
    ! --------------------------------------------------
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction_char = unique_directions%directions_char(j)
      if (direction_char=='x') then
        displacement = [1.0_dp,0.0_dp,0.0_dp]
      elseif (direction_char=='y') then
        displacement = [0.0_dp,1.0_dp,0.0_dp]
      elseif (direction_char=='z') then
        displacement = [0.0_dp,0.0_dp,1.0_dp]
      else
        call err()
      endif
      
      atom_string = left_pad(atom, str(maxval(unique_directions%atoms)))
      paths = [ sdir//'/atom.'//atom_string//'.+d'//direction_char, &
              & sdir//'/atom.'//atom_string//'.-d'//direction_char  ]
      
      ! Make harmonic run directories.
      do k=1,2
        call mkdir(paths(k))
        
        ! Move relevant atom.
        if(k==1) then
          call supercell%atoms(atom)%set_cartesian_position( &
             & supercell%atoms(atom)%cartesian_position() + displacement)
        elseif(k==2) then
          call supercell%atoms(atom)%set_cartesian_position( &
             & supercell%atoms(atom)%cartesian_position() - 2*displacement)
        endif
          
        ! Write dft input files.
        input_filename = make_input_filename(file_type,seedname)
        call StructureData_to_input_file(     &
                   & file_type,               &
                   & supercell,               &
                   & wd//'/'//input_filename, &
                   & paths(k)//'/'//input_filename)
      enddo
      
      ! Reset moved atom.
      call supercell%atoms(atom)%set_cartesian_position( &
         & supercell%atoms(atom)%cartesian_position() + displacement)
    enddo
  enddo
end subroutine
end module
