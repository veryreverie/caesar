! ======================================================================
! The second stage of Caesar.
! Runs DFT for harmonic calculations.
! ======================================================================
module run_harmonic_module
  use constants_module, only : dp
  use string_module
  use io_module
contains

! ----------------------------------------------------------------------
! Generate keywords and helptext.
! ----------------------------------------------------------------------
function run_harmonic_keywords() result(keywords)
  use keyword_module
  implicit none
  
  type(KeywordData), allocatable :: keywords(:)
  
  keywords = [                                                                &
  & KeywordData( 'supercells_to_run',                                        &
  &               'supercells_to_run is the first and last supercell to run. &
  &These should be specified as two integers separated by spaces.'),          &
  & KeywordData( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1'),                                         &
  & KeywordData( 'run_script',                                               &
  &               'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.) ]
end function

function run_harmonic_mode() result(output)
  use caesar_modes_module
  implicit none
  
  type(CaesarMode) :: output
  
  output%mode_name = 'run_harmonic'
  output%description = 'Runs DFT calculations set up by setup_harmonic. &
     &should be run after setup_harmonic.'
  output%keywords = run_harmonic_keywords()
  output%main_subroutine => run_harmonic
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_harmonic(arguments)
  use setup_harmonic_module
  use unique_directions_module
  use dictionary_module
  use IFile_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory
  type(String) :: wd
  
  ! Previous user inputs.
  type(Dictionary) :: setup_harmonic_arguments
  type(IFile)      :: no_supercells_file
  integer          :: no_supercells
  type(String)     :: dft_code
  type(String)     :: seedname
  
  ! Terminal inputs.
  integer      :: supercells_to_run(2)
  integer      :: no_cores
  type(String) :: run_script
  
  ! Atom and direction data.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  type(String)           :: atom_string
  
  ! Temporary variables.
  integer      :: i,j,k
  character(1) :: signs(2)
  type(String) :: dir
  type(String) :: sdir
  integer      :: result_code
  
  ! --------------------------------------------------
  ! Get inputs from user.
  ! --------------------------------------------------
  wd = arguments%value('working_directory')
  supercells_to_run = int(split(arguments%value('supercells_to_run')))
  no_cores = int(arguments%value('no_cores'))
  run_script = arguments%value('run_script')
  
  ! --------------------------------------------------
  ! Read in arguments to previous calculations.
  ! --------------------------------------------------
  no_supercells_file = wd//'/no_supercells.dat'
  no_supercells = int(no_supercells_file%line(1))
  
  setup_harmonic_arguments = setup_harmonic_keywords()
  call setup_harmonic_arguments%read_file(wd//'/setup_harmonic.used_settings')
  dft_code = setup_harmonic_arguments%value('dft_code')
  seedname = setup_harmonic_arguments%value('seedname')
  
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
  ! Run calculations
  ! --------------------------------------------------
  do i=supercells_to_run(1),supercells_to_run(2)
    sdir = wd//'/Supercell_'//left_pad(i,str(no_supercells))
    
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction = unique_directions%directions_char(j)
      atom_string = left_pad(atom,str(maxval(unique_directions%atoms)))
      
      signs = [ '+', '-' ]
      do k=1,2
        dir = sdir//'/atom.'//atom_string//'.'//signs(k)//'d'//direction
        call print_line('')
        call print_line('Running calculation in directory '//dir)
        result_code = system_call( 'cd '//wd//'; '//run_script//' '// &
           & dft_code//' '//dir//' '//no_cores//' '//seedname)
        call print_line('Result code: '//result_code)
      enddo
    enddo
  enddo
end subroutine
end module
