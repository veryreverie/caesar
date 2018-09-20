! ======================================================================
! The second stage of Caesar.
! Runs DFT for harmonic calculations.
! ======================================================================
module run_harmonic_module
  use common_module
  
  use setup_harmonic_module
  implicit none
  
  private
  
  public :: run_harmonic
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_harmonic() result(output)
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'run_harmonic'
  output%description = 'Runs DFT calculations set up by setup_harmonic. &
     &should be run after setup_harmonic.'
  output%keywords = [                                                         &
     & KeywordData( 'supercells_to_run',                                      &
     &              'supercells_to_run is the first and last supercell to &
     &run. These should be specified as two integers separated by spaces.'),  &
     & KeywordData( 'run_script',                                             &
     &              'run_script is the path to the script for running DFT. An &
     &example run script can be found in doc/input_files.',                   &
     &              is_path=.true.),                                          &
     & KeywordData( 'no_cores',                                               &
     &              'no_cores is the number of cores on which DFT will be &
     &run. This is passed to the specified run script.',                      &
     &              default_value='1'),                                       &
     & KeywordData( 'calculation_type',                                       &
     &              'calculation_type specifies whether any electronic &
     &structure calculations should be run in addition to the user-defined &
     &script. Settings are: "none" and "quip".',                              &
     &              default_value='none') ]
  output%main_subroutine => run_harmonic_subroutine
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_harmonic_subroutine(arguments)
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! User inputs.
  integer      :: supercells_to_run(2)
  integer      :: no_cores
  type(String) :: run_script
  type(String) :: calculation_type
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  integer          :: no_supercells
  type(String)     :: file_type
  type(String)     :: seedname
  
  ! Electronic structure calculation runner.
  type(CalculationRunner) :: calculation_runner
  
  ! Atom and direction data.
  type(UniqueDirection), allocatable :: unique_directions(:)
  integer                            :: atom
  type(String)                       :: direction
  type(String)                       :: atom_string
  
  ! Files and Directories.
  type(IFile)  :: structure_file
  type(IFile)  :: no_supercells_file
  type(IFile)  :: unique_directions_file
  type(String) :: supercell_dir
  type(String) :: run_dir
  
  ! Temporary variables.
  integer :: i,j,ialloc
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  supercells_to_run = int(split_line(arguments%value('supercells_to_run')))
  run_script = arguments%value('run_script')
  no_cores = int(arguments%value('no_cores'))
  calculation_type = arguments%value('calculation_type')
  
  ! --------------------------------------------------
  ! Read in arguments to previous calculations.
  ! --------------------------------------------------
  setup_harmonic_arguments = Dictionary(setup_harmonic())
  call setup_harmonic_arguments%read_file('setup_harmonic.used_settings')
  file_type = setup_harmonic_arguments%value('file_type')
  seedname = setup_harmonic_arguments%value('seedname')
  
  ! --------------------------------------------------
  ! Read in previously calculated data.
  ! --------------------------------------------------
  no_supercells_file = IFile('no_supercells.dat')
  no_supercells = int(no_supercells_file%line(1))
  
  ! --------------------------------------------------
  ! Check user inputs.
  ! --------------------------------------------------
  if (supercells_to_run(1)<=0) then
    call print_line('')
    call print_line('Error: first supercell must be > 0')
    call err()
  elseif (supercells_to_run(1)>no_supercells) then
    call print_line('')
    call print_line('Error: first supercell must be <= '//no_supercells)
    call err()
  elseif (supercells_to_run(2)<supercells_to_run(1)) then
    call print_line('')
    call print_line('Error: first supercell must be <= last supercell.')
    call err()
  elseif (supercells_to_run(2)>no_supercells) then
    call print_line('')
    call print_line('Error: last supercell must be <= '//no_supercells)
    call err()
  elseif (no_cores<=0) then
    call print_line('')
    call print_line('Error: no. cores must be >= 1')
    call err()
  elseif (.not. file_exists(run_script)) then
    call print_line('')
    call print_line('Error: '//run_script//' does not exist.')
    call err()
  endif
  
  ! --------------------------------------------------
  ! Initialise calculation runner.
  ! --------------------------------------------------
  calculation_runner = CalculationRunner(  &
    & file_type         = file_type,       &
    & seedname          = seedname,        &
    & run_script        = run_script,      &
    & no_cores          = no_cores,        &
    & calculation_type  = calculation_type )
  
  ! --------------------------------------------------
  ! Run calculations
  ! --------------------------------------------------
  ! Loop over supercells.
  do i=supercells_to_run(1),supercells_to_run(2)
    supercell_dir = 'Supercell_'//left_pad(i,str(no_supercells))
    
    unique_directions_file = IFile(supercell_dir//'/unique_directions.dat')
    unique_directions = UniqueDirection(unique_directions_file%sections())
    
    ! Loop over displacements within each supercell.
    do j=1,size(unique_directions)
      atom = unique_directions(j)%atom_id
      direction = unique_directions(j)%direction
      atom_string = left_pad(atom, str(maxval(unique_directions%atom_id)))
      
      run_dir = supercell_dir//'/atom.'//atom_string//'.'//direction
      
      ! Run calculation at each displacement.
      call calculation_runner%run_calculation(run_dir)
    enddo
    
    deallocate(unique_directions, stat=ialloc); call err(ialloc)
  enddo
end subroutine
end module
