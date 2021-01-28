! ======================================================================
! The first stage of Caesar.
! Generates supercells, and prepares harmonic DFT calculations.
! ======================================================================
module caesar_setup_harmonic_module
  use caesar_common_module
  
  use caesar_harmonic_data_module
  use caesar_generate_supercells_module
  implicit none
  
  private
  
  public :: startup_setup_harmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
subroutine startup_setup_harmonic()
  implicit none
  
  type(CaesarMode) :: mode
  
  mode%mode_name = 'setup_harmonic'
  mode%description = 'Sets up harmonic calculation. Generates supercells, &
     &and prepares DFT inputs.'
  mode%keywords = [                                                           &
     & KeywordData( 'file_type',                                              &
     &              'file_type is the file type which will be used for &
     &single-point energy calculations. Settings are: "castep", &
     &"quantum_espresso", "caesar" and "xyz".',                               &
     &              default_value='castep'),                                  &
     & KeywordData( 'seedname',                                               &
     &              'seedname is the seedname from which file names are &
     &constructed.'),                                                         &
     & KeywordData( 'q-point_grid',                                           &
     &              'q-point_grid is the number of q-points in each direction &
     &in a Monkhorst-Pack grid. This should be specified as three integers &
     &separated by spaces.'),                                                 &
     & KeywordData( 'symmetry_precision',                                     &
     &              'symmetry_precision is the tolerance up to which a&
     &symmetry is accepted. In order for a symmetry to be accepted, it must &
     &transform the position of every atom to within symmetry_precision of an &
     &atom of the same element. There must be no other atom within &
     &symmetry_precision of this point. symmetry_precision should be given in &
     &Bohr. Ideally, symmetry_precision should be much smaller than the &
     &minimum inter-atomic distance, and much larger than the geometry &
     &optimisation tolerance.',                                               &
     &              default_value='0.01'),                                    &
     & KeywordData( 'harmonic_displacement',                                  &
     &              'harmonic_displacement is the distance in bohr by which &
     &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
     &surface.',                                                              &
     &              default_value='0.01'),                                    &
     & KeywordData( 'snap_to_symmetry',                                       &
     &              'snap_to_symmetry specifies whether or not to enforce &
     &exact symmetry on the structure before running calculations. If &
     &snap_to_symmetry is set to false, and the symmetry error is large, this &
     &may lead to errors and potential crashes in later calculations.',       &
     &              default_value='true'),                                    &
     & KeywordData( 'loto_direction',                                         &
     &              'loto_direction specifies the direction (in reciprocal &
     &co-ordinates from which the gamma point is approached when calculating &
     &LO/TO corrections. loto_direction should be specified as three &
     &fractions, e.g. 1/2 1/2 1/2. This does not have to be normalised. If &
     &loto_direction is not set, then LO/TO corrections are not calculated. &
     &N.B. LO/TO corrections can currently only be calculated from Castep &
     &DFTB calculations, and only if the q-point grid is 1x1x1. If &
     &loto_direction is specified here, it will be used for the whole &
     &calculation and may not be re-specified in calculate_normal_modes or &
     &calculate_potential.',                                                  &
     &              is_optional = .true.)                                     ]
  mode%main_subroutine => setup_harmonic_subroutine
  
  call add_mode(mode)
end subroutine

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine setup_harmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User input variables.
  type(String)         :: file_type
  type(String)         :: seedname
  integer              :: grid(3)
  real(dp)             :: symmetry_precision
  real(dp)             :: harmonic_displacement
  logical              :: snap_to_symmetry
  type(Fractionvector) :: loto_direction
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
  
  ! Common harmonic data.
  type(HarmonicData)            :: harmonic_data
  type(StructureData)           :: structure
  type(StructureData)           :: large_supercell
  type(QpointData), allocatable :: qpoints(:)
  
  ! Non-diagonal supercell data.
  integer                          :: no_supercells
  type(StructureData), allocatable :: supercells(:)
  
  ! Directories.
  type(String) :: supercell_dir
  type(String) :: calculation_dir
  
  ! Perturbation direction information.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(CartesianDisplacement)        :: displacement
  type(StructureData)                :: displaced_structure
  
  ! Files.
  type(OFile) :: structure_file
  type(Ofile) :: large_supercell_file
  type(OFile) :: no_supercells_file
  type(OFile) :: qpoints_file
  type(OFile) :: supercell_file
  type(OFile) :: unique_directions_file
  type(OFile) :: calculation_directories_file
  
  ! Temporary variables.
  integer :: i,j
  
  ! Read in arguments from user.
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  grid = int(split_line(arguments%value('q-point_grid')))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  snap_to_symmetry = lgcl(arguments%value('snap_to_symmetry'))
  if (arguments%is_set('loto_direction')) then
    loto_direction = FractionVector(arguments%value('loto_direction'))
    if (any(grid/=1)) then
      call print_line(ERROR//': Unable to calculate LO/TO splitting if &
         &q-point grid is not 1x1x1.')
      call quit()
    endif
  endif
  
  ! Initialise calculation writer.
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  ! Calculate harmonic data. This includes:
  !    - The primitive cell structure.
  !    - The large supercell structure.
  !    - The q-points.
  if (arguments%is_set('loto_direction')) then
    harmonic_data = HarmonicData( file_type,          &
                                & seedname,           &
                                & grid,               &
                                & symmetry_precision, &
                                & snap_to_symmetry,   &
                                & loto_direction      )
  else
    harmonic_data = HarmonicData( file_type,          &
                                & seedname,           &
                                & grid,               &
                                & symmetry_precision, &
                                & snap_to_symmetry    )
  endif
  
  structure = harmonic_data%structure
  large_supercell = harmonic_data%large_supercell
  qpoints = harmonic_data%qpoints
  
  ! Generate non-diagonal supercells.
  if (arguments%is_set('loto_direction')) then
    supercells = generate_supercells( structure,          &
                                    & qpoints,            &
                                    & symmetry_precision, &
                                    & loto_direction      )
  else
    supercells = generate_supercells(structure, qpoints, symmetry_precision)
  endif
  no_supercells = size(supercells)
  
  ! Write out structure, q-point and supercell data.
  structure_file = OFile('structure.dat')
  call structure_file%print_lines(structure)
  
  large_supercell_file = OFile('large_supercell.dat')
  call large_supercell_file%print_lines(large_supercell)
  
  qpoints_file = OFile('qpoints.dat')
  call qpoints_file%print_lines(qpoints,separating_line='')
  
  no_supercells_file = OFile('no_supercells.dat')
  call no_supercells_file%print_line(no_supercells)
  
  ! Loop over supercells, writing out calculation directories for each.
  do i=1,no_supercells
    supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
    
    call mkdir(supercell_dir)
    
    ! Write out structure files.
    supercell_file = OFile(supercell_dir//'/structure.dat')
    call supercell_file%print_lines(supercells(i))
    
    ! Calculate which forces need calculating.
    unique_directions = calculate_unique_directions( supercells(i),        &
                                                   & harmonic_displacement )
    unique_directions_file = OFile(              &
       & supercell_dir//'/unique_directions.dat' )
    call unique_directions_file%print_lines( unique_directions, &
                                           & separating_line='')
    
    ! Write energy and force calculation input files.
    do j=1,size(unique_directions)
      ! Construct displaced structure.
      displacement = CartesianDisplacement(unique_directions(j), supercells(i))
      displaced_structure = displace_structure(supercells(i), displacement)
      
      ! Write calculation input files.
      calculation_dir = supercell_dir//'/atom.'                            // &
                      & left_pad( unique_directions(j)%atom_id,               &
                      &           str(maxval(unique_directions%atom_id)) ) // &
                      & '.'//unique_directions(j)%direction
      call calculation_writer%write_calculation( displaced_structure, &
                                               & calculation_dir      )
    enddo
  enddo
  
  ! Write out calculation directories to file.
  calculation_directories_file = OFile('harmonic_calculation_directories.dat')
  call calculation_directories_file%print_lines( &
     & calculation_writer%directories_written()  )
end subroutine
end module
