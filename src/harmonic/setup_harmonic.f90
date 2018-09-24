! ======================================================================
! The first stage of Caesar.
! Generates supercells, and prepares harmonic DFT calculations.
! ======================================================================
module setup_harmonic_module
  use common_module
  
  use generate_supercells_module
  implicit none
  
  private
  
  public :: setup_harmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function setup_harmonic() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'setup_harmonic'
  output%description = 'Sets up harmonic calculation. Generates supercells, &
     &and prepares DFT inputs.'
  output%keywords = [                                                         &
  & KeywordData( 'file_type',                                                 &
  &              'file_type is the file type which will be used for &
  &single-point energy calculations. Settings are: "castep", "caesar" and &
  &"xyz". Support for xyz files requires Caesar to be linked against Quip.',  &
  &              default_value='castep'),                                     &
  & KeywordData( 'seedname',                                                  &
  &              'seedname is the seedname from which file names are &
  &constructed.'),                                                            &
  & KeywordData( 'q-point_grid',                                              &
  &              'q-point_grid is the number of q-points in each direction &
  &in a Monkhorst-Pack grid. This should be specified as three integers &
  &separated by spaces.'),                                                    &
  & KeywordData( 'symmetry_precision',                                        &
  &              'In order for a symmetry to be accepted, it must transform &
  &the position of every atom to within symmetry_precision of an atom of the &
  &same element. symmetry_precision should be given in Bohr.',                &
  &              default_value='0.1'),                                        &
  & KeywordData( 'harmonic_displacement',                                     &
  &              'harmonic_displacement is the distance in bohr by which &
  &atoms will be displaced when mapping the harmonic Born-Oppenheimer &
  &surface.',                                                                 &
  &              default_value='0.01') ]
  output%main_subroutine => setup_harmonic_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine setup_harmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User input variables.
  type(String) :: file_type
  type(String) :: seedname
  
  ! Electronic structure calculation writer.
  type(CalculationWriter) :: calculation_writer
  
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
  type(String) :: supercell_dir
  type(String) :: calculation_dir
  
  ! Perturbation direction information.
  type(UniqueDirection), allocatable :: unique_directions(:)
  type(CartesianDisplacement)        :: displacement
  type(StructureData)                :: displaced_structure
  
  ! Files.
  type(String) :: input_filename
  type(OFile)  :: structure_file
  type(Ofile)  :: large_supercell_file
  type(OFile)  :: no_supercells_file
  type(OFile)  :: qpoints_file
  type(OFile)  :: supercell_file
  type(OFile)  :: unique_directions_file
  
  ! Temporary variables.
  integer :: i,j
  
  ! --------------------------------------------------
  ! Get settings from user, and check them.
  ! --------------------------------------------------
  file_type = arguments%value('file_type')
  seedname = arguments%value('seedname')
  grid = int(split_line(arguments%value('q-point_grid')))
  symmetry_precision = dble(arguments%value('symmetry_precision'))
  harmonic_displacement = dble(arguments%value('harmonic_displacement'))
  
  ! --------------------------------------------------
  ! Initialise calculation writer.
  ! --------------------------------------------------
  calculation_writer = CalculationWriter( file_type = file_type, &
                                        & seedname  = seedname   )
  
  ! --------------------------------------------------
  ! Read in input files.
  ! --------------------------------------------------
  input_filename = make_input_filename(file_type, seedname)
  structure = input_file_to_StructureData(file_type, input_filename)
  
  ! --------------------------------------------------
  ! Generate symmetries of structure.
  ! --------------------------------------------------
  if (size(structure%symmetries)==0) then
    call structure%calculate_symmetry(symmetry_precision)
  endif
  
  ! --------------------------------------------------
  ! Generate large supercell, for which all q-points are G-vectors.
  ! --------------------------------------------------
  large_supercell_matrix = mat( [ grid(1), 0      , 0     ,    &
                              &   0      , grid(2), 0     ,    &
                              &   0      , 0      , grid(3) ], &
                              & 3,3)
  large_supercell = construct_supercell( structure,             &
                                       & large_supercell_matrix )
  
  ! --------------------------------------------------
  ! Generate q-points in Monkhorst-Pack grid.
  ! --------------------------------------------------
  qpoints = generate_qpoints(large_supercell)
  
  ! --------------------------------------------------
  ! Generate non-diagonal supercells.
  ! --------------------------------------------------
  supercells = generate_supercells(structure, qpoints, symmetry_precision)
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
  
  ! Loop over supercells.
  do i=1,no_supercells
    supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
    
    call mkdir(supercell_dir)
    
    ! Write out structure files.
    supercell = supercells(i)
    supercell_file = OFile(supercell_dir//'/structure.dat')
    call supercell_file%print_lines(supercell)
    
    ! Calculate which forces need calculating.
    unique_directions = calculate_unique_directions( supercell,            &
                                                   & harmonic_displacement )
    unique_directions_file = OFile(              &
       & supercell_dir//'/unique_directions.dat' )
    call unique_directions_file%print_lines( unique_directions, &
                                           & separating_line='')
    
    ! --------------------------------------------------
    ! Write energy and force calculation input files.
    ! --------------------------------------------------
    do j=1,size(unique_directions)
      ! Construct displaced structure.
      displacement = unique_directions(j)%cartesian_displacement(supercell)
      displaced_structure = displace_structure(supercell,displacement)
      
      ! Write calculation input files.
      calculation_dir = supercell_dir//'/atom.'                            // &
                      & left_pad( unique_directions(j)%atom_id,               &
                      &           str(maxval(unique_directions%atom_id)) ) // &
                      & '.'//unique_directions(j)%direction
      call calculation_writer%write_calculation( displaced_structure, &
                                               & calculation_dir      )
    enddo
  enddo
end subroutine
end module
