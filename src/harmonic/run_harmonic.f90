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
  use help_module
  implicit none
  
  type(KeywordData) :: keywords(3)
  
  keywords = [                                                                &
  & make_keyword( 'supercells_to_run',                                        &
  &               'supercells_to_run is the first and last supercell to run. &
  &These should be specified as two integers separated by spaces.'),          &
  & make_keyword( 'no_cores',                                                 &
  &               'no_cores is the number of cores on which DFT will be run. &
  &This is passed to the specified run script.',                              &
  &               default_value='1'),                                         &
  & make_keyword( 'run_script',                                               &
  &               'run_script is the path to the script for running DFT. An &
  &example run script can be found in doc/input_files.',                      &
  &               is_path=.true.) ]
end function

! ----------------------------------------------------------------------
! Main program.
! ----------------------------------------------------------------------
subroutine run_harmonic(arguments)
  use unique_directions_module
  use dictionary_module
  implicit none
  
  type(Dictionary), intent(in) :: arguments
  
  ! Working directory
  type(String) :: wd
  
  ! Previous user inputs.
  type(Dictionary)          :: setup_harmonic_arguments
  type(String), allocatable :: no_supercells_file(:)
  integer                   :: no_supercells
  type(String)              :: dft_code
  type(String)              :: seedname
  
  ! Terminal inputs.
  integer      :: supercells_to_run(2)
  integer      :: no_cores
  type(String) :: run_script
  
  ! Atom and direction data.
  type(UniqueDirections) :: unique_directions
  integer                :: atom
  character(1)           :: direction
  
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
  no_supercells_file = read_lines(wd//'/no_supercells.dat')
  no_supercells = int(no_supercells_file(1))
  
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
    sdir = wd//'/Supercell_'//i
    
    unique_directions = read_unique_directions_file( &
       & sdir//'/unique_directions.dat')
    
    do j=1,size(unique_directions)
      atom = unique_directions%atoms(j)
      direction = unique_directions%directions_char(j)
      
      signs = [ '+', '-' ]
      do k=1,2
        dir = sdir//'/atom.'//atom//'.'//signs(k)//'d'//direction
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
